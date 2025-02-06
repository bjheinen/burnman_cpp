/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_UTIL_EOS_HPP_INCLUDED
#define BURNMAN_UTIL_EOS_HPP_INCLUDED

/**
 * Sruct to hold excesses data for EOS modifiers
 */
struct Excesses {
    double G = 0.0;
    double dGdT = 0.0;
    double dGdP = 0.0;
    double d2GdT2 = 0.0;
    double d2GdP2 = 0.0;
    double d2GdPdT = 0.0;
    // Overloaded += operator to add Excesses
    Excesses& operator+=(const Excesses& other) {
        G += other.G;
        dGdT += other.dGdT;
        dGdP += other.dGdP;
        d2GdT2 += other.d2GdT2;
        d2GdP2 += other.d2GdP2;
        d2GdPdT += other.d2GdPdT;
        return *this;
    }
};

#endif // BURNMAN_UTIL_EOS_HPP_INCLUDED