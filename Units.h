#ifndef UNITS_H_INCLUDED
#define UNITS_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <random>
#include <time.h>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <assert.h>


#define kB 8.617330e-5
#define e 1.6e-19
#define pi 3.14159265359

#define repetition 8

enum Direction {x, y, z};

enum AtomType { C, O, S, H};

enum FragmentType { edot, tos};

AtomType mass2Type(float masse);

struct Bound
{
    int m_id1;
    int m_id2;
};

struct Vect3
{
    double m_x;
    double m_y;
    double m_z;

    Vect3(): m_x(0), m_y(0), m_z(0)
    {
    }

	Vect3(double x, double y, double z):m_x(x), m_y(y), m_z(z)
	{
	}

    Vect3& operator=(const Vect3& v)
    {
		m_x=v.m_x;
		m_y=v.m_y;
		m_z=v.m_z;
		return *this;
	}

	Vect3& operator+=(const Vect3& v)
    {
		m_x+=v.m_x;
		m_y+=v.m_y;
		m_z+=v.m_z;
		return *this;
	}

    Vect3& operator-=(const Vect3& v)
    {
		m_x-=v.m_x;
		m_y-=v.m_y;
		m_z-=v.m_z;
		return *this;
	}

    Vect3 operator/(double r) const
    {
		return Vect3(m_x/r, m_y/r, m_z/r);
	}

    double norm() const
    {
		return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
	}

	double operator*(const Vect3 &r) const
	{
		return (m_x*r.m_x + m_y*r.m_y + m_z*r.m_z);
	}

    Vect3 operator-(const Vect3 &r) const
    {
		return Vect3(m_x-r.m_x, m_y-r.m_y, m_z-r.m_z);
	}

    Vect3 operator+(const Vect3 &r) const
    {
		return Vect3(m_x+r.m_x, m_y+r.m_y, m_z+r.m_z);
	}

	Vect3 operator*(double r) const
	{
		return Vect3(m_x*r,m_y*r,m_z*r);
	}

	void toZero()
	{
        m_x = 0;
        m_y = 0;
        m_z = 0;
	}

	bool isNull()
	{
	    return (m_x==0 && m_y==0 && m_z==0);
	}

	// Cross product
	Vect3 operator^(const Vect3 &r) const
	{
		return Vect3( m_y * r.m_z - m_z * r.m_y, m_z * r.m_x - m_x * r.m_z, m_x * r.m_y - m_y * r.m_x );
	}

	static double distance(const Vect3 &v1, const Vect3 &v2)
	{
        return ( sqrt( (v1.m_x-v2.m_x)*(v1.m_x-v2.m_x)  + (v1.m_y-v2.m_y)*(v1.m_y-v2.m_y) + (v1.m_z-v2.m_z)*(v1.m_z-v2.m_z) ) );
	}

	// Distance with periodic boundary conditions
    static double distancePbc(const Vect3 &v1, const Vect3 &v2, const Vect3 &boxSize, bool &periodic)
	{
	    periodic = false;
	    double dxabs(fabs(v2.m_x-v1.m_x));
	    double dx = (std::min(dxabs, boxSize.m_x - dxabs ));
	    if (dx == boxSize.m_x - dxabs) {periodic = true;}
	    double dyabs(fabs(v2.m_y-v1.m_y));
	    double dy = (std::min(dyabs, boxSize.m_y - dyabs ));
	    if (dy == boxSize.m_y - dyabs) {periodic = true;}
        double dzabs(fabs(v2.m_z-v1.m_z));
	    double dz = (std::min(dzabs, boxSize.m_z - dzabs ));
	    if (dz == boxSize.m_z - dzabs) {periodic = true;}
        return ( sqrt( dx*dx + dy * dy + dz *dz ) );
	}

    // Displacement vector with periodic boundary conditions
    static Vect3 displacementPbc(const Vect3 &v1, const Vect3 &v2, const Vect3 &boxSize)
	{
	    Vect3 v;

	    double dxabs(fabs(v2.m_x-v1.m_x));
	    if(dxabs < boxSize.m_x - dxabs)
        {
            v.m_x = v2.m_x - v1.m_x;
        }
        else
        {
            if ( v2.m_x > v1.m_x ) { v.m_x = -(boxSize.m_x - dxabs); }
            else { v.m_x = (boxSize.m_x - dxabs); }

        }

	    double dyabs(fabs(v2.m_y-v1.m_y));
        if(dyabs < boxSize.m_y - dyabs)
        {
            v.m_y = v2.m_y - v1.m_y;
        }
        else
        {
            if ( v2.m_y > v1.m_y) { v.m_y = -(boxSize.m_y - dyabs); }
            else { v.m_y = (boxSize.m_y - dyabs); }

        }

        double dzabs(fabs(v2.m_z-v1.m_z));
        if(dzabs < boxSize.m_z - dzabs)
        {
            v.m_z = v2.m_z - v1.m_z;
        }
        else
        {
            if ( v2.m_z > v1.m_z ) { v.m_z = -(boxSize.m_z - dzabs); }
            else { v.m_z = (boxSize.m_z - dzabs); }

        }
        return v;
	}


	// Displacement vector with periodic boundary conditions + indication of the concerned pbc (x,y, or/and z in Vect3 pbc, with 0 or 1 or -1 depending on the direction)
    static Vect3 displacementPbc(const Vect3 &v1, const Vect3 &v2, const Vect3 &boxSize, Vect3 &pbc)
	{
	    Vect3 v;

	    double dxabs(fabs(v2.m_x-v1.m_x));
	    if(dxabs < boxSize.m_x - dxabs)
        {
            v.m_x = v2.m_x - v1.m_x;
            pbc.m_x = 0;
        }
        else
        {
            if ( v2.m_x > v1.m_x ) { v.m_x = -(boxSize.m_x - dxabs); pbc.m_x = -1; }
            else { v.m_x = (boxSize.m_x - dxabs); pbc.m_x = 1; }

        }

	    double dyabs(fabs(v2.m_y-v1.m_y));
        if(dyabs < boxSize.m_y - dyabs)
        {
            v.m_y = v2.m_y - v1.m_y;
            pbc.m_y = 0;
        }
        else
        {
            if ( v2.m_y > v1.m_y) { v.m_y = -(boxSize.m_y - dyabs); pbc.m_y = -1;}
            else { v.m_y = (boxSize.m_y - dyabs); pbc.m_y = 1;}

        }

        double dzabs(fabs(v2.m_z-v1.m_z));
        if(dzabs < boxSize.m_z - dzabs)
        {
            v.m_z = v2.m_z - v1.m_z;
            pbc.m_z = 0;
        }
        else
        {
            if ( v2.m_z > v1.m_z ) { v.m_z = -(boxSize.m_z - dzabs); pbc.m_z = -1;}
            else { v.m_z = (boxSize.m_z - dzabs); pbc.m_z = 1;}

        }
        return v;
	}


		// Displacement vector between two sites, with periodic boundary conditions, given the indication of the concerned pbc (x,y, or/and z in Vect3 pbc, with 0 or 1 or -1 depending on the direction)
    static Vect3 siteDisplacementPbc(const Vect3 &v1, const Vect3 &v2, const Vect3 &boxSize, const Vect3 &pbc)
	{
	    Vect3 v;

	    double dxabs(fabs(v2.m_x-v1.m_x));
	    if(pbc.m_x == 0)
        {
            v.m_x = v2.m_x - v1.m_x;
        }
        else
        {
            if ( pbc.m_x == -1 ) { v.m_x = -(boxSize.m_x - dxabs); }
            else { v.m_x = (boxSize.m_x - dxabs); }

        }

	    double dyabs(fabs(v2.m_y-v1.m_y));
        if(pbc.m_y == 0)
        {
            v.m_y = v2.m_y - v1.m_y;
        }
        else
        {
            if ( pbc.m_y == -1) { v.m_y = -(boxSize.m_y - dyabs); }
            else { v.m_y = (boxSize.m_y - dyabs); }

        }

        double dzabs(fabs(v2.m_z-v1.m_z));
        if(pbc.m_z == 0)
        {
            v.m_z = v2.m_z - v1.m_z;
        }
        else
        {
            if ( pbc.m_z == -1 ) { v.m_z = -(boxSize.m_z - dzabs);}
            else { v.m_z = (boxSize.m_z - dzabs);}

        }
        return v;
	}
};


bool operator==(Vect3 const& v1, Vect3 const& v2);

struct Atom
{
    int m_name;
    AtomType m_type;
    Vect3 m_pos;
};

struct RigidFragment
{
    FragmentType m_type;
    int m_molNbr; // To what polymer chain it belongs
    int m_siteNbr; // To what conjugated segment it belongs
    int m_fragNbr; // Index of the fragment in the conjugated segment

    Vect3 m_center;
    Vect3 m_orientN;
    Vect3 m_orientP;
    Vect3 m_orientT;

    std::vector<Atom> m_atom;

    void computeGeom(bool invert);

};

struct Polymer
{
    int m_chainLength;
    std::vector<RigidFragment> fragment;
};


struct Topo
{
    int m_chainLength;
    int m_pedotNbr;
    int m_tosNbr;
    std::vector<int> m_HTerm;
    std::vector<AtomType> m_atomPedot;
    std::vector<Bound> m_boundPedot;
    std::vector<std::vector<int> > m_neighborPedot;
    std::vector<AtomType> m_atomTos;
    std::vector<Bound> m_boundTos;
    std::vector<std::vector<int> > m_neighborTos;
    std::vector<Bound> m_topoMap;

    void readTopo(std::string fileName);
    void createTopoMap();
};

bool crossingRaySphere(const Vect3 &rf1,const Vect3 &rf2,const Vect3 &rfInfluence,const double &sphereR, const Vect3 &boxSize);


#endif // UNITS_H_INCLUDED
