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
// It offers the subset of the igzstream interface the line-based loaders
// use (open, fail, clear, close, getline), so it is a drop-in replacement
// there.

#ifndef __GZLINEREADER_H__
#define __GZLINEREADER_H__

#include <zlib.h>
#include <cstring>
#include <string>
#include <vector>

class GzLineReader {
public:
    GzLineReader() : buf(1 << 20) {}
    explicit GzLineReader(const char *filename) : GzLineReader() { open(filename); }
    ~GzLineReader() { close(); }
    GzLineReader(const GzLineReader &) = delete;
    GzLineReader &operator=(const GzLineReader &) = delete;

    void open(const char *filename) {
        close();
        f = gzopen(filename, "rb");
        bad = (f == nullptr);
        if (f) gzbuffer(f, 1 << 17);
        pos = len = 0;
        atEof = false;
    }
    bool fail() const { return bad; }
    bool is_open() const { return f != nullptr; }
    void clear() { bad = false; }
    void close() {
        if (f) { gzclose(f); f = nullptr; }
        pos = len = 0;
        atEof = false;
    }

    bool getline(std::string &out) {
        out.clear();
        bool extracted = false;
        for (;;) {
            if (pos < len) {
                const char *s = buf.data() + pos;
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
            if (atEof || f == nullptr) {
                if (!extracted) bad = true;   // like a stream at end of file
                return extracted;
            }
            int n = gzread(f, buf.data(), static_cast<unsigned>(buf.size()));
            if (n <= 0) {
                atEof = true;
                continue;
            }
            pos = 0;
            len = static_cast<size_t>(n);
        }
    }

private:
    gzFile f = nullptr;
    std::vector<char> buf;
    size_t pos = 0, len = 0;
    bool bad = false;
    bool atEof = false;
};

inline bool getline(GzLineReader &r, std::string &line) { return r.getline(line); }

#endif
