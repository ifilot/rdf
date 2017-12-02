#ifndef _POSITION_H
#define _POSITION_H

class Position {
public:
    float x,y,z;

    Position():
    x(0), y(0), z(0) {}

    Position(float _x, float _y, float _z):
    x(_x), y(_y), z(_z) {}
};

#endif //_POSITION_H
