//
//  numbers2strings.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 23/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "numbers2strings.hpp"

/*
 * C++ version std::string style "itoa":
 //http://www.jb.man.ac.uk/~slowe/cpp/itoa.html
 */
std::string itoa(int value, unsigned int base) {
    const char digitMap[] = "0123456789abcdef";
    std::string buf;
    // Guard:
    if (base == 0 || base > 16) {
        // Error: may add more trace/log output here
        return buf;
    }
    // Take care of negative int:
    std::string sign;
    int _value = value;
    // Check for case when input is zero:
    if (_value == 0) return "0";
    if (value < 0) {
        _value = -value;
        sign = "-";
    }
    
    // Translating number to string with base:
    for (int i = 30; _value && i ; --i) {
        buf = digitMap[ _value % base ] + buf;
        _value /= base;
    }
    return sign.append(buf);
}
std::string double2string(double dbl){
    std::ostringstream strs;
    strs << dbl;
    return strs.str();
}
