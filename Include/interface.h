#ifndef INTERFACE_H
#define INTERFACE_H

#include <string>

enum class object_type {PathwaySet, D2Grid, Distribution};

class Interface
{
public:
    Interface();
    object_type ObjectType;


};

#endif // INTERFACE_H
