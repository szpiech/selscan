// gzlinereader.h -- buffered line reader for gzip-compressed or plain text.
//
// P1: igzstream reads through a 303-byte stream buffer and std::getline
// extracts one character at a time through it, which made VCF loading
// about 20x slower than decompression. GzLineReader reads 1 MiB blocks
// with gzread (zlib reads uncompressed files transparently) and splits
// lines with memchr. getline() has the semantics of std::getline: the
// '\n' is removed, any '\r' is kept, a last line without '\n' is
// returned, and false means nothing was extracted.
//
// M9: a second thread runs gzread ahead of the parser, filling up to
// three 1 MiB blocks, so decompression overlaps the parsing of the lines.
// The calling thread only splits lines. Block order is preserved, so the
// lines are exactly those of the sequential reader. The thread is detached
// and shares the file and blocks through a shared_ptr, so an error that
// calls exit() part way through a file leaves no thread to join.
//
// It offers the subset of the igzstream interface the line-based loaders
// use (open, fail, clear, close, getline), so it is a drop-in replacement
// there.

#ifndef __GZLINEREADER_H__
#define __GZLINEREADER_H__

#include <zlib.h>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

class GzLineReader {
public:
    GzLineReader() {}
    explicit GzLineReader(const char *filename) : GzLineReader() { open(filename); }
    ~GzLineReader() { close(); }
    GzLineReader(const GzLineReader &) = delete;
    GzLineReader &operator=(const GzLineReader &) = delete;

    void open(const char *filename) {
        close();
        gzFile f = gzopen(filename, "rb");
        bad = (f == nullptr);
        if (!f) return;
        gzbuffer(f, 1 << 17);
        st = std::make_shared<State>();
        st->f = f;
        for (int i = 0; i < NBLOCKS; i++) { st->blocks[i].assign(BLOCK, 0); st->freeList.push_back(i); }
        std::thread(produce, st).detach();
    }
    bool fail() const { return bad; }
    bool is_open() const { return st != nullptr; }
    void clear() { bad = false; }
    void close() {
        if (st) {
            { std::lock_guard<std::mutex> lk(st->m); st->stop = true; }
            st->cv.notify_all();
            st.reset();   // the last owner (here or the thread) closes the file
        }
        pos = len = 0;
        cur = -1;
        atEof = false;
    }

    bool getline(std::string &out) {
        out.clear();
        bool extracted = false;
        for (;;) {
            if (pos < len) {
                const char *s = st->blocks[cur].data() + pos;
                const char *nl = static_cast<const char *>(std::memchr(s, '\n', len - pos));
                if (nl) {
                    out.append(s, nl - s);
                    pos += (nl - s) + 1;
                    return true;
                }
                out.append(s, len - pos);
                pos = len;
                extracted = true;
            }
            if (atEof || !st) {
                if (!extracted) bad = true;   // like a stream at end of file
                return extracted;
            }
            if (!nextBlock()) {
                atEof = true;
                continue;
            }
        }
    }

private:
    static constexpr int NBLOCKS = 3;
    static constexpr size_t BLOCK = 1 << 20;

    // shared by the reader and its thread
    struct State {
        gzFile f = nullptr;
        std::vector<char> blocks[NBLOCKS];
        std::deque<std::pair<int, size_t>> filled;   // block index, bytes
        std::deque<int> freeList;
        std::mutex m;
        std::condition_variable cv;
        bool stop = false, producerDone = false;
        ~State() { if (f) gzclose(f); }
    };

    // thread: decompress into free blocks, in order, until end of file or stop
    static void produce(std::shared_ptr<State> s) {
        for (;;) {
            int b;
            {
                std::unique_lock<std::mutex> lk(s->m);
                s->cv.wait(lk, [&] { return s->stop || !s->freeList.empty(); });
                if (s->stop) return;
                b = s->freeList.front();
                s->freeList.pop_front();
            }
            int n = gzread(s->f, s->blocks[b].data(), static_cast<unsigned>(BLOCK));
            {
                std::lock_guard<std::mutex> lk(s->m);
                if (n > 0) s->filled.push_back({b, static_cast<size_t>(n)});
                else { s->freeList.push_back(b); s->producerDone = true; }
            }
            s->cv.notify_all();
            if (n <= 0) return;
        }
    }

    // reader: hand the current block back and take the next filled one
    bool nextBlock() {
        std::unique_lock<std::mutex> lk(st->m);
        if (cur >= 0) { st->freeList.push_back(cur); cur = -1; st->cv.notify_all(); }
        st->cv.wait(lk, [&] { return !st->filled.empty() || st->producerDone; });
        if (st->filled.empty()) { pos = len = 0; return false; }
        cur = st->filled.front().first;
        len = st->filled.front().second;
        st->filled.pop_front();
        pos = 0;
        return true;
    }

    std::shared_ptr<State> st;
    int cur = -1;
    size_t pos = 0, len = 0;
    bool bad = false;
    bool atEof = false;
};

inline bool getline(GzLineReader &r, std::string &line) { return r.getline(line); }

#endif
