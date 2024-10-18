#ifndef NUMERICS_UTILITY_HPP
#define NUMERICS_UTILITY_HPP


float f(float x) {
    return x * x - 2.0f;
}

float bisection(float (*f)(float x), float l, float r) {
    float fl = f(l); 
    float x = 0.0f; float fx = 0.0f;
    int iter_max = 10;  // use tolerance also?
    for (int i = 0; i < iter_max; ++i) {
        x = (l + r) / 2.0f;
        fx = f(x); 
        cout << "x and f(x): " << x << ", " << fx << endl;;
        cout << "  cur [ " << l << ", " << r << "]" << endl;
        if ( (fl > 0 && fx < 0) || ( fl < 0 && fx > 0 ) ) {
           r = x;
        } else {
           l = x;
        }
        cout << "  new [ " << l << ", " << r << "]" << endl;
    }
    std::cout << "final answer: " << x << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    return x;
}

float false_position(float (*f)(float x), float l, float r) {
    float fl = f(l); float fr = f(r);
    float x = 0.0f; float fx = 0.0f;
    int iter_max = 10;  // use tolerance also?
    for (int i = 0; i < iter_max; ++i) { 
        float m = (fr - fl)/(r - l);    
        float b = fl - m * l;
        x = - b / m;
        fx = f(x);
        cout << "x and f(x): " << x << ", " << fx << endl;;
        cout << "  cur [ " << l << ", " << r << "]" << endl;
        if ( (fl > 0 && fx < 0) || ( fl < 0 && fx > 0 ) ) {
           r = x;
        } else {
           l = x;
        }
        cout << "  new [ " << l << ", " << r << "]" << endl;
    }
    std::cout << "final answer: " << x << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    return x;
}

#endif
