
#include <iostream>
#include <string>
#include <cmath>
#include<iomanip>
#include "../hapmap/bitset.h"
using namespace std;
// Macro to define a test case
#define TEST_CASE(test_name) void test_name()

// Macro to check a condition and report success or failure
#define CHECK(condition) \
    do { \
        if (!(condition)) { \
            std::cerr << "TEST FAILED: " << #condition << " in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            return; \
        } else {\
            std::cout << "TEST PASSED: " << #condition << " in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        } \ 
    } while (0)

// Function to run a test case
void run_test(void (*test)(), const std::string& test_name) {
    std::cout << "Running TEST " << test_name << "..." << std::endl;
    test();
    std::cout  <<"Test"<< " ended: " << test_name  << std::endl;
}

TEST_CASE(test_bitset) {
    MyBitset *bp = new MyBitset(10);
    MyBitset b = *bp;

    MyBitset *bp2 = new MyBitset(10);
    MyBitset b2= *bp2;
    b.set_bit(3);
    b.set_bit(9);
    b.set_bit(10);

    for(int locus_after_filter = 0; locus_after_filter < b.nbits; locus_after_filter++){
        if(locus_after_filter==0){
            MyBitset* b1 = &b;
            for (int k = 0; k < b1->nwords; k++) {
                ->bits[k] = b1->bits[k] ;
            }
            hapEntries[locus_after_filter].xorbitset->num_1s = b1->num_1s;
        }else{
            MyBitset* b1 =(hapEntries[locus_after_filter].hapbitset);
            MyBitset* b2 = (hapEntries[locus_after_filter-1].hapbitset);

            int sum = 0;
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[locus_after_filter].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                sum += __builtin_popcountll(hapEntries[locus_after_filter].xorbitset->bits[k]);
            }
            hapEntries[locus_after_filter].xorbitset->num_1s = sum;
            
        }
    }

    for(int i = 0; i<10; i++){
        if(i==3 ||  i==9 || i==10)
            CHECK(b.get_bit(i));
        else
            CHECK(!b.get_bit(i));
    }

    vector<int> p = b.get_bits();
    for(int i = 0; i<p.size(); i++){
       cout<<p[i]<<endl;
    }
}

TEST_CASE(test_addition) {
    int a = 2;
    int b = 3;
    int result = a + b;
    CHECK(result == 5);
}


int main (int argc, char* argv[]){

    run_test(test_addition, "test_addition");
    run_test(test_bitset, "test_bitset");


    return 0;
}

