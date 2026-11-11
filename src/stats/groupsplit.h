// groupsplit.h -- reusable group-split bookkeeping for the EHH kernels.
//
// P5: at every step of every core the kernels collected, for each current
// haplotype group, the haplotypes carrying one allele at the new site in an
// unordered_map<int, vector<int>> (old group id -> member haplotypes), then
// split the groups. That is one hash node and one vector allocation per
// group touched per step; on nSL with 2000 haplotypes it was most of the
// run time.
//
// GroupSplitMap offers the subset of that interface the kernels use --
// m[g].push_back(h), m.clear(), range-for over (first = g, second =
// members), second.size(), second[0], range-for over second -- with flat
// storage that is reused: no allocation once the buffers have grown.
//
// Differences from unordered_map: groups are visited in the order they
// were first touched (unordered_map's order was unspecified), members of a
// group keep insertion order as before. The kernels' updates are sums of
// pair counts, which are integers held exactly in doubles, and new group
// ids are only labels, so the visiting order does not change results.
//
// threadSplitMap(slot) returns a cleared per-thread instance (slot 0, and
// slot 1 for the XP pooled map), so the buffers are also reused across
// cores. A kernel must not hold a slot across a call to another kernel.

#ifndef __GROUPSPLIT_H__
#define __GROUPSPLIT_H__

#include <cstddef>
#include <vector>

class GroupSplitMap {
public:
    class Members {   // members of one group, contiguous
    public:
        Members(const int *b, const int *e) : b_(b), e_(e) {}
        size_t size() const { return static_cast<size_t>(e_ - b_); }
        const int *begin() const { return b_; }
        const int *end() const { return e_; }
        int operator[](size_t k) const { return b_[k]; }
    private:
        const int *b_, *e_;
    };
    struct Entry { int first; Members second; };

    class Adder {
    public:
        Adder(GroupSplitMap &m, int g) : m_(m), g_(g) {}
        void push_back(int h) { m_.add(g_, h); }
    private:
        GroupSplitMap &m_;
        int g_;
    };

    class const_iterator {
    public:
        const_iterator(const GroupSplitMap *m, size_t k) : m_(m), k_(k) {}
        Entry operator*() const {
            const int *base = m_->members.data();
            return Entry{m_->touched[k_], Members(base + m_->groupBegin[k_], base + m_->groupBegin[k_ + 1])};
        }
        const_iterator &operator++() { ++k_; return *this; }
        bool operator!=(const const_iterator &o) const { return k_ != o.k_; }
    private:
        const GroupSplitMap *m_;
        size_t k_;
    };

    Adder operator[](int g) { return Adder(*this, g); }

    void add(int g, int h) {
        if (static_cast<size_t>(g) >= cnt.size()) { cnt.resize(g + 1, 0); cursor.resize(g + 1); }
        if (cnt[g]++ == 0) touched.push_back(g);
        pairs.push_back(g);
        pairs.push_back(h);
        built = false;
    }
    bool empty() const { return touched.empty(); }
    size_t size() const { return touched.size(); }
    void clear() {
        for (int g : touched) cnt[g] = 0;
        touched.clear();
        pairs.clear();
        built = false;
    }

    const_iterator begin() const { build(); return const_iterator(this, 0); }
    const_iterator end() const { return const_iterator(this, touched.size()); }

private:
    // Counting sort of the (group, member) pairs into contiguous runs,
    // groups in first-touched order, members in insertion order.
    void build() const {
        if (built) return;
        groupBegin.resize(touched.size() + 1);
        size_t off = 0;
        for (size_t k = 0; k < touched.size(); k++) {
            groupBegin[k] = off;
            cursor[touched[k]] = off;
            off += cnt[touched[k]];
        }
        groupBegin[touched.size()] = off;
        members.resize(off);
        for (size_t j = 0; j < pairs.size(); j += 2) members[cursor[pairs[j]]++] = pairs[j + 1];
        built = true;
    }

    std::vector<int> cnt;        // members per group id this step (0 when untouched)
    std::vector<int> touched;    // group ids in first-touched order
    std::vector<int> pairs;      // group, member, group, member, ...
    mutable std::vector<size_t> cursor, groupBegin;
    mutable std::vector<int> members;
    mutable bool built = false;
};

inline GroupSplitMap &threadSplitMap(int slot) {
    static thread_local GroupSplitMap maps[2];
    GroupSplitMap &m = maps[slot];
    m.clear();
    return m;
}

#endif
