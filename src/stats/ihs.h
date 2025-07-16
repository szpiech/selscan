#ifndef __SELSCAN_IHS_H__
#define __SELSCAN_IHS_H__

#include "selscan-stats.h"
#include "../thread_pool.h"
#include <unordered_map>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

#define ACTION_ON_ALL_SET_BITS(hapbitset, ACTION)         \
    for (int k = 0; k < (hapbitset->nwords); k++) {             \
        uint64_t bitset = (hapbitset->bits)[k];                 \
        while (bitset != 0) {                        \
            uint64_t t = bitset & -bitset;           \
            int r = __builtin_ctzl(bitset);          \
            int set_bit_pos = (k * 64 + r);          \
            bitset ^= t;                             \
            ACTION;                                  \
        }                                            \
    }

struct IhhComponents {
    double derived_right;
    double derived_left;
    double ancestral_right;
    double ancestral_left;
};

struct OutputUnphasedIHH {
    double iHH2;
    double iHH0;
    double ciHH2;
    double ciHH0;
};


class IHS: public SelscanStats{
    public:
        IHS(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){  

        }
        void main(); 
        pair<double, double> calc_ihh1(int locus);  
        IhhComponents calc_ihh1_details(int locus);
        pair<double, double> infer_missing(int locus);  
        std::mutex mutex_log;
        //std::unique_ptr<std::mutex> mutex_log = std::make_unique<std::mutex>();

    protected:
        void updateEHH_from_split_unphased( unordered_map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, bool* is1, bool* is2, int* group_core);
        string getOrder(uint64_t n_c2, uint64_t n_c1, uint64_t n_c0);


        /// assume already validated ranges
        std::vector<std::pair<int, int>> parse_ranges(const std::string& range_str) {
            std::vector<std::pair<int, int>> ranges;
            std::stringstream ss(range_str);
            std::string token;

            while (std::getline(ss, token, ',')) {
                size_t dash_pos = token.find('-');
                // if (dash_pos == std::string::npos) {
                // HANDLE_ERROR_NOLOG("Invalid range format: " + token);
                // }

                int start = std::stoi(token.substr(0, dash_pos));
                int end = std::stoi(token.substr(dash_pos + 1));

                    //if (start > end) std::swap(start, end);
                ranges.emplace_back(start, end);

                // try {
                //     int start = std::stoi(token.substr(0, dash_pos));
                //     int end = std::stoi(token.substr(dash_pos + 1));

                //     if (start > end) std::swap(start, end);
                //     ranges.emplace_back(start, end);
                // } catch (const std::invalid_argument&) {
                //     HANDLE_ERROR_NOLOG("Invalid number in range: " + token);
                // } catch (const std::out_of_range&) {
                //     HANDLE_ERROR_NOLOG("Number out of range in range: " + token);
                // }
            }

            return ranges;
        }

        void create_directories(const std::string& path) {
            std::istringstream iss(path);
            std::string token;
            std::string current_path;

            while (std::getline(iss, token, '/')) {
                current_path += token + "/";
                mkdir(current_path.c_str(), 0755);  // ignores errors if exists
            }
            
        }

        std::string basename(const std::string& path) {
            size_t pos = path.find_last_of("/\\");
            if (pos == std::string::npos) return path;
            return path.substr(pos + 1);
        }
        bool is_position_in_ranges(int pos, const std::vector<pair<int, int> >& ranges) {
            for (const auto& r : ranges) {
                if (pos >= r.first && pos <= r.second) {
                    return true;
                }
            }
            return false;
        }

    private:
        //static pthread_mutex_t mutex_log;
        int max_extend;

        //phased_ihs
        pair<double, double> calc_ehh_unidirection(int locus, bool downstream);

        //unphased_ihs  
        OutputUnphasedIHH calc_ehh_unidirection_unphased(int locus, bool downstream);
};

#endif
