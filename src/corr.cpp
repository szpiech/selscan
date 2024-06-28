#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// Template function to write the contents of a vector to a file
template <typename T>
void writeXYToFile(const std::vector<T>& x, const std::vector<T>& y, const std::string& filename) {
    std::ofstream file(filename);  // Open the file in write mode
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (int i = 0; i< x.size(); i++) {
        const T& elem_x = x[i];
        const T& elem_y = y[i];

        file << elem_x << " "<< elem_y << std::endl;  // Write each element followed by a newline
    }

    file.close();  // Close the file
}



template <typename T>
std::vector<T> readColumn(const std::string& filename, int column, bool skipHeader) {
     std::vector<T> column_values;
    std::ifstream file(filename);  // Open the file
    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        return column_values;
    }

   
    std::string line;

    if(skipHeader) {
        std::getline(file, line);
    }
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
        int column_index = 0;
        T column_value;

        while (iss >> value) {
            column_index++;
            if (column_index == column) {
                try {
                    column_value = std::stod(value);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid value encountered: " << value << std::endl;
                    continue;
                } catch (const std::out_of_range& e) {
                    std::cerr << "Value out of range: " << value << std::endl;
                    continue;
                }
                column_values.push_back(column_value);
                break;
            }
        }
    }


    file.close();
    return column_values;
}


double mean(const std::vector<double>& v) {
    double sum = 0;
    for (double val : v) {
        sum += val;
    }
    return sum / v.size();
}

double covariance(const std::vector<double>& x, const std::vector<double>& y) {
    double mean_x = mean(x);
    double mean_y = mean(y);
    double cov = 0;

    for (size_t i = 0; i < x.size(); ++i) {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return cov / x.size();
}

double variance(const std::vector<double>& v) {
    double mean_v = mean(v);
    double var = 0;

    for (double val : v) {
        var += (val - mean_v) * (val - mean_v);
    }
    return var / v.size();
}


template <typename T>
std::vector<T> findIntersection(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> intersection;
    
    // Copy the input vectors to avoid modifying the original vectors
    std::vector<T> sorted_v1 = v1;
    std::vector<T> sorted_v2 = v2;

    // Sort the vectors // assume sorted
    
    // Use two pointers to find the intersection
    size_t i = 0, j = 0;
    while (i < sorted_v1.size() && j < sorted_v2.size()) {
        if (sorted_v1[i] < sorted_v2[j]) {
            i++;
        } else if (sorted_v1[i] > sorted_v2[j]) {
            j++;
        } else {
            intersection.push_back(sorted_v1[i]);
            i++;
            j++;
        }
    }
    
    return intersection;
}



double pearson_correlation(const std::vector<double>& x, const std::vector<double>& y) {
    double cov = covariance(x, y);
    double std_x = std::sqrt(variance(x));
    double std_y = std::sqrt(variance(y));
    return cov / (std_x * std_y);
}

int main() {
    std::vector<double> hb_ihh1;
    std::vector<double> sb_ihh1;

    std::vector<double> hb_ihh0;
    std::vector<double> sb_ihh0;

    std::vector<int> hb_id;
    std::vector<int> sb_id;

    hb_id = readColumn<int>("out/hapbin_chr10.txt", 1, true);
    sb_id = readColumn<int>("out/selbin_chr10.txt", 1, false);

    hb_ihh1 = readColumn<double>("out/hapbin_chr10.txt", 5, true);
    sb_ihh1 = readColumn<double>("out/selbin_chr10.txt", 5, false);

    hb_ihh0 = readColumn<double>("out/hapbin_chr10.txt", 4, true);
    sb_ihh0 = readColumn<double>("out/selbin_chr10.txt", 6, false);

    // hb_ihh0 = readColumn<double>("out/hapbin_chr10.txt", 6, true);
    // sb_ihh0 = readColumn<double>("out/selbin_chr10.txt", 7, false);
    //HAPBIN 6 - IHS
    //selbin 7 - IHS

    std::vector<int> common_ids = findIntersection(hb_id, sb_id);
    std::vector<double> hb_ihh1_intersect;
    std::vector<double> sb_ihh1_intersect;
    hb_ihh1_intersect.reserve(common_ids.size());
    sb_ihh1_intersect.reserve(common_ids.size());
    
    std::vector<double> hb_ihh0_intersect;
    std::vector<double> sb_ihh0_intersect;
    hb_ihh0_intersect.reserve(common_ids.size());
    sb_ihh0_intersect.reserve(common_ids.size());


    int it_hb = 0;
    int it_sb = 0;
    for(int id : common_ids) {
        while(hb_id[it_hb] < id) {
            it_hb++;
        }
        if(it_hb >= hb_id.size()) {
            break;
        }
        if(hb_id[it_hb] > id){
            continue;
        }
        while(sb_id[it_sb] < id) {
            it_sb++;
        }
        if(it_sb >= sb_id.size()) {
            break;
        }
        if(sb_id[it_sb] > id){
            continue;
        }
        if(sb_id[it_sb] == hb_id[it_hb] ){
            hb_ihh1_intersect.push_back(hb_ihh1[it_hb]);
            sb_ihh1_intersect.push_back(sb_ihh1[it_sb]);

            hb_ihh0_intersect.push_back(hb_ihh0[it_hb]);
            sb_ihh0_intersect.push_back(sb_ihh0[it_sb]);
        }
       

        //std::cout << id << " ";
    }

    int cnt = 0;
    for (double val : hb_ihh1_intersect) {
        std::cout << val << " ";
        cnt++;
        if(cnt > 20){
            break;
        }
    }
    std::cout << std::endl;

    cnt = 0;
    for (double val : sb_ihh1_intersect) {
        std::cout << val << " ";
        cnt++;
        if(cnt > 20){
            break;
        }
    }
    std::cout << std::endl;


    std::vector<double>& x = hb_ihh0_intersect;
    std::vector<double>& y = sb_ihh0_intersect;
    if (x.size() != y.size()) {
        std::cerr << "Error: Datasets must be of the same size." << std::endl;
        return 1;
    }

    double correlation = pearson_correlation(x, y);
    writeXYToFile(hb_ihh0_intersect, sb_ihh0_intersect, "out/ihh0_hb_sb.txt");
    std::cout << "Pearson correlation coefficient: " << correlation << std::endl;
    //6th column of selbin. 4 5 ihh1 ihh0
    //hapbin 4 5 ihh0 ihh1
    //selbin

    //index by col 1 both

    
    return 0;
}
