// normio.h -- buffered token input and formatted output for selscan norm.
//
// P6: the read and normalise passes used operator>> for every field and
// ostream insertion for every number written; on a 5M-row iHS file that
// was ~8 of ~11 s (locale lookups, stream sentries and per-character
// string growth). These classes keep the same semantics:
//
// NormTokenReader returns whitespace-delimited tokens exactly as
// operator>> into a std::string would (tokens may continue on the next
// line, empty lines are skipped, an exhausted file yields an empty token),
// reading whole lines through GzLineReader. num() applies the B18 rule of
// SelscanNorm::parseNum: a token that is not entirely a finite number
// gives NaN.
//
// NormOut appends to a 1 MiB buffer. putG() formats a double with
// "%.6g", the conversion libc++ and libstdc++ use for ostream's default
// floatfield and precision 6, so output is byte-identical to the
// previous `fout << x`. close() reports write errors, which the ofstream
// code ignored.

#ifndef __NORMIO_H__
#define __NORMIO_H__

#include "gzlinereader.h"
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

class NormTokenReader {
public:
    explicit NormTokenReader(const std::string &filename) : r(filename.c_str()) {}
    bool fail() const { return r.fail(); }

    // Whole line, discarding anything left of the current line.
    bool getline(std::string &out) { p = line.size(); return r.getline(out); }

    // Number of whitespace-delimited fields on the next line; that line
    // becomes the source of the following tokens.
    int fieldsInNextLine() {
        if (!r.getline(line)) line.clear();
        p = 0;
        int nf = 0;
        bool in = false;
        for (char c : line) {
            bool sp = std::isspace(static_cast<unsigned char>(c)) != 0;
            if (!sp && !in) nf++;
            in = !sp;
        }
        return nf;
    }

    bool next(const char *&s, size_t &n) {
        for (;;) {
            const size_t L = line.size();
            while (p < L && std::isspace(static_cast<unsigned char>(line[p]))) p++;
            if (p < L) {
                size_t b = p;
                while (p < L && !std::isspace(static_cast<unsigned char>(line[p]))) p++;
                s = line.data() + b;
                n = p - b;
                return true;
            }
            if (!r.getline(line)) { line.clear(); p = 0; s = line.data(); n = 0; return false; }
            p = 0;
        }
    }
    void str(std::string &out) { const char *s; size_t n; next(s, n); out.assign(s, n); }
    void skip(int k) { const char *s; size_t n; for (int i = 0; i < k; i++) next(s, n); }
    long integer() {   // as operator>> into an int for a well-formed integer token
        const char *s; size_t n;
        if (!next(s, n)) return 0;
        return std::strtol(s, nullptr, 10);
    }
    double num() {
        const char *s; size_t n;
        next(s, n);
        if (n == 0) return NAN;
        char *e = nullptr;
        double v = std::strtod(s, &e);   // the token is followed by whitespace or '\0'
        if (e != s + n || !std::isfinite(v)) return NAN;
        return v;
    }

private:
    GzLineReader r;
    std::string line;
    size_t p = 0;
};

class NormOut {
public:
    NormOut() { b.reserve(BUF + 256); }
    ~NormOut() { if (f) { flush(); std::fclose(f); } }
    NormOut(const NormOut &) = delete;
    NormOut &operator=(const NormOut &) = delete;

    bool open(const std::string &filename) { f = std::fopen(filename.c_str(), "w"); return f != nullptr; }
    void put(const char *s, size_t n) { b.append(s, n); if (b.size() >= BUF) flush(); }
    void put(const std::string &s) { put(s.data(), s.size()); }
    void put(char c) { b.push_back(c); if (b.size() >= BUF) flush(); }
    void putInt(long v) { char t[24]; int k = std::snprintf(t, sizeof t, "%ld", v); put(t, k); }
    void putG(double v) { char t[40]; int k = std::snprintf(t, sizeof t, "%.6g", v); put(t, k); }
    // false if any write or the close failed
    bool close() {
        if (!f) return !err;
        flush();
        if (std::fclose(f) != 0) err = true;
        f = nullptr;
        return !err;
    }

private:
    static const size_t BUF = 1 << 20;
    void flush() {
        if (f && !b.empty() && std::fwrite(b.data(), 1, b.size(), f) != b.size()) err = true;
        b.clear();
    }
    FILE *f = nullptr;
    std::string b;
    bool err = false;
};

#endif
