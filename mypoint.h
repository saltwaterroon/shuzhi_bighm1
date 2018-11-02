//
//  Header.h
//  MyFirstOpenCV
//
//  Created by h-xiao16 on 2018/10/30.
//  Copyright © 2018年 h-xiao16. All rights reserved.
//


#define EPSILON 0.00001f
#include <cmath>


// =========================================
// 3-mypoint
// =========================================
class mypoint
{
public:
    
    // Position
    float x, y, pixel;
    
    // Default constructor
    mypoint()
    : x( 0 ), y( 0 ), pixel( 0 ) {}
    
    // Element constructor
    mypoint( float x, float y, float pixel )
    : x( x ), y( y ), pixel( pixel ) {}
    
    // Copy constructor
    mypoint( const mypoint& a )
  : x( a.x ), y( a.y ), pixel( a.pixel ) {}
    
    // Norm (len^2)
    inline float len() const { return x*x + y*y; }
    
    // distance between tow points
    inline float distance() const { return (float)sqrt(len()); }
    
    mypoint &operator += ( const mypoint &src ) { x += src.x; y += src.y; pixel += src.pixel; return *this; }
    mypoint operator + ( const mypoint &src ) const { mypoint tmp( *this ); return ( tmp += src ); }
    mypoint &operator -= ( const mypoint &src ) { x -= src.x; y -= src.y; pixel -= src.pixel; return *this; }
    mypoint operator - ( const mypoint &src ) const { mypoint tmp( *this ); return ( tmp -= src ); }
    
    mypoint operator - () const { return mypoint(-x,-y,-pixel); }
    
    mypoint &operator *= ( const float src ) { x *= src; y *= src; pixel *= src;  return *this; }
    mypoint operator * ( const float src ) const { mypoint tmp( *this ); return ( tmp *= src ); }
    mypoint &operator /= ( const float src ) { x /= src; y /= src; pixel /= src; return *this; }
    mypoint operator / ( const float src ) const { mypoint tmp( *this ); return ( tmp /= src ); }
    
    bool operator == ( const mypoint& b) const { return ((*this)-b).len() < EPSILON; }
    //bool operator == ( const mypoint& b) const { return x==b.x && y==b.y && pixel==b.pixel; }
};

// Left hand float multplication
inline mypoint operator * ( const float src, const mypoint& v ) { mypoint tmp( v ); return ( tmp *= src ); }

// Dot product
inline float dot( const mypoint& a, const mypoint& b )
{ return a.x*b.x + a.y*b.y; }


