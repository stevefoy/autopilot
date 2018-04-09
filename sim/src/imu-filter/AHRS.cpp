/* -*- indent-tabs-mode:T; c-basic-offset:8; tab-width:8; -*- vi: set ts=8:
 *  $Id: AHRS.cpp,v 2.16 2003/03/15 05:53:38 tramm Exp $
 *
 * (c) Trammell Hudson
 * (c) Aaron Kahn
 *
 * AHRS simulator based on Kalman filtering of the gyro and
 * accelerometer data.  Converted from Aaron's matlab code
 * to use the C++ math library.
 *
 **************
 *
 *  This file is part of the autopilot simulation package.
 *
 *  Autopilot is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Autopilot is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Autopilot; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <AHRS.h>

#include <mat/Vector.h>
#include <mat/Matrix.h>
#include <mat/Matrix_Invert.h>
#include <mat/Quat.h>
#include <mat/Nav.h>
#include <mat/Conversions.h>
#include <mat/Kalman.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>


#include "macros.h"

namespace imufilter
{

using namespace util;
using namespace libmat;
using namespace std;


void
AHRS::make_a_matrix(
	Matrix<N,N> &		A,
	const Vector<3> &	pqr
) const
{
	const state_t &		x = this->state;
	const double		p = pqr[0] / 2.0;
	const double		q = pqr[1] / 2.0;
	const double		r = pqr[2] / 2.0;

	const double		q0 = x.q[0] / 2.0;
	const double		q1 = x.q[1] / 2.0;
	const double		q2 = x.q[2] / 2.0;
	const double		q3 = x.q[3] / 2.0;

	A[0][0] =  0;
	A[0][1] = -p;
	A[0][2] = -q;
	A[0][3] = -r;
	A[0][4] =  q1;
	A[0][5] =  q2;
	A[0][6] =  q3;

	A[1][0] =  p;
	A[1][1] =  0;
	A[1][2] =  r;
	A[1][3] = -q;
	A[1][4] = -q0;
	A[1][5] =  q3;
	A[1][6] = -q2;

	A[2][0] =  q;
	A[2][1] = -r;
	A[2][2] =  0;
	A[2][3] =  p;
	A[2][4] = -q3;
	A[2][5] = -q0;
	A[2][6] =  q1;

	A[3][0] =  r;
	A[3][1] =  q;
	A[3][2] = -p;
	A[3][3] =  0;
	A[3][4] =  q2;
	A[3][5] = -q1;
	A[3][6] = -q0;

	for( int i=4 ; i<N ; i++ )
		for( int j=0 ; j<N ; j++ )
			A[i][j] = 0;
}


void
AHRS::propagate_state(
	const Vector<3> &	pqr
)
{
	state_t &		x = this->state;
	const double		p = pqr[0] / 2.0;
	const double		q = pqr[1] / 2.0;
	const double		r = pqr[2] / 2.0;

	const double		q0 = x.q[0];
	const double		q1 = x.q[1];
	const double		q2 = x.q[2];
	const double		q3 = x.q[3];

	Vector<4>		Xdot(
		-p * q1 - q * q2 - r * q3,
		 p * q0 - q * q3 + r * q2,
		 p * q3 + q * q0 - r * q1,
		-p * q2 + q * q1 + r * q0
	);

	Xdot *= this->dt;

	x.q += Xdot;
	x.q.norm_self();
}


void
AHRS::propagate_covariance(
	const Matrix<N,N> &	A
)
{
	Matrix<N,N>		Pdot( this->Q );
	Pdot += A * this->P;
	Pdot += this->P * A.transpose();
	Pdot *= this->dt;

	this->P += Pdot;

	this->trace = 0;

	for( int i=0 ; i<N ; i++ )
		this->trace += this->P[i][i];	
}



template<
	int			m
>
void
AHRS::do_kalman(
	const Matrix<m,N> &	C,
	const Matrix<m,m> &	R,
	const Vector<m> &	eTHETA
)
{
	state_t &		x( this->state );

	// We throw away the K result
	Matrix<N,m>		K;

	// Kalman() wants a vector, not an object.  Serialize the
	// state data into this vector, then extract it out again
	// once we're done with the loop.
	Vector<N>		X_vect;

	X_vect[0]	= x.q[0];
	X_vect[1]	= x.q[1];
	X_vect[2]	= x.q[2];
	X_vect[3]	= x.q[3];

	X_vect[4]	= x.bias[0];
	X_vect[5]	= x.bias[1];
	X_vect[6]	= x.bias[2];

	Kalman(
		this->P,
		X_vect,
		C,
		R,
		eTHETA,
		K
	);

	x.q[0]		= X_vect[0];
	x.q[1]		= X_vect[1];
	x.q[2]		= X_vect[2];
	x.q[3]		= X_vect[3];

	x.bias[0]	= X_vect[4];
	x.bias[1]	= X_vect[5];
	x.bias[2]	= X_vect[6];

	x.q.norm_self();

}


void
AHRS::kalman_attitude_update(
	const Vector<3> &	accel,
	const Matrix<3,3> &	DCM,
	const Vector<3> &	THETAe
)
{
	double			err;
	state_t &		x( this->state );

	// compute the euler angles from the accelerometers
	const Vector<3>		THETAm( accel2euler( accel, THETAe[2] ) );

	// make the C matrix
	Matrix<2,N>		C;

	// PHI section
	err = 2.0 / ( sqr(DCM[2][2]) + sqr(DCM[1][2]) );

	C[0][0] = ( x.q[1] * DCM[2][2]                            ) * err;
	C[0][1] = ( x.q[0] * DCM[2][2] + 2.0 * x.q[1] * DCM[1][2] ) * err;
	C[0][2] = ( x.q[3] * DCM[2][2] + 2.0 * x.q[2] * DCM[1][2] ) * err;
	C[0][3] = ( x.q[2] * DCM[2][2]                            ) * err;
	C[0][4] = 0.0;
	C[0][5] = 0.0;
	C[0][6] = 0.0;

	// THETA section
	err = -1.0 / sqrt(1.0 - sqr(DCM[0][2]) );

	C[1][0] = -2.0 * x.q[2] * err;
	C[1][1] =  2.0 * x.q[3] * err;
	C[1][2] = -2.0 * x.q[0] * err;
	C[1][3] =  2.0 * x.q[1] * err;
	C[1][4] = 0.0;
	C[1][5] = 0.0;
	C[1][6] = 0.0;


	// compute the error; this should be ( THETAm - THETAe ),
	// but we can only use the pitch and roll angles here
	Vector<2>		eTHETA;
	eTHETA[0] = THETAm[0] - THETAe[0];
	eTHETA[1] = THETAm[1] - THETAe[1];

	this->do_kalman(
		C,
		this->R_attitude,
		eTHETA
	);
}


void
AHRS::kalman_compass_update(
	double			heading,		// degrees
	const Matrix<3,3> &	DCM,
	const Vector<3> &	THETAe
)
{
	state_t &		x( this->state );

	const double		q0 = x.q[0];
	const double		q1 = x.q[1];
	const double		q2 = x.q[2];
	const double		q3 = x.q[3];

	Matrix<1,N>		C( 0 );

	// PSI section
	const double		err = 2 / (sqr(DCM[0][0]) + sqr(DCM[0][1]));
	C[0][0] = err * ( q3 * DCM[0][0] );
	C[0][1] = err * ( q2 * DCM[0][0] );
	C[0][2] = err * ( q1 * DCM[0][0] + 2.0 * q2 * DCM[0][1] );
	C[0][3] = err * ( q0 * DCM[0][0] + 2.0 * q3 * DCM[0][1] );

	// Compute the error
	Vector<1>		eTHETA;

	eTHETA[0] = heading - THETAe[2];
	if( eTHETA[0] > C_PI )
		eTHETA[0] -= 2.0 * C_PI;
	else
	if( eTHETA[0] < -C_PI )
		eTHETA[0] += 2.0 * C_PI;

	this->do_kalman(
		C,
		this->R_heading,
		eTHETA
	);
}


AHRS::AHRS(
	double			dt
) :
	dt( dt )
{
	this->reset();
}


void
AHRS::reset()
{
	//this->P			= eye<N,double>();
	this->P.fill();
	P[4][4] = 1;
	P[5][5] = 1;
	P[6][6] = 1;

	this->Q.fill();
	this->R_attitude.fill();
	this->R_heading.fill();

	// Quaterion attitude estimate noise
	this->Q[0][0] = 0.0001;
	this->Q[1][1] = 0.0001;
	this->Q[2][2] = 0.0001;
	this->Q[3][3] = 0.0001;

	// Gyro bias
	this->Q[4][4] = 0.03;
	this->Q[5][5] = 0.03;
	this->Q[6][6] = 0.03;

	this->R_attitude[0][0] = 0.3;
	this->R_attitude[1][1] = 0.3;

	this->R_heading[0][0] = 0.5;
}


/*
 * We assume that the vehicle is still during the first sample
 * and use the values to help us determine the zero point for the
 * gyro bias and accelerometers.
 *
 * You must call this once you have the samples from the IMU
 * and compass.  Perhaps throw away the first few to let things
 * stabilize.
 */
void
AHRS::initialize(
	const Vector<3> &	accel,
	const Vector<3> &	pqr,
	double			heading
)
{
	this->state.bias	= pqr;
	this->state.q		= euler2quat( accel2euler( accel, heading ) );
	this->theta		= quat2euler( this->state.q );
}


void
AHRS::imu_update(
	const Vector<3> &	accel,
	const Vector<3> &	pqr_raw
)
{
	const state_t &		x( this->state );
	const Vector<3>		pqr( pqr_raw - x.bias );
	Matrix<N,N>		A;

	this->make_a_matrix( A, pqr );
	this->propagate_state( pqr );
	this->propagate_covariance( A );

	/* Compute the DCM and angle estimate here */
	const Matrix<3,3> 	DCM( quatDC( x.q ) );
	const Vector<3>		THETAe( quat2euler( x.q ) );

	this->kalman_attitude_update(
		accel,
		DCM,
		THETAe
	);

	/* compute angles from quaternions */
	this->theta	= quat2euler( x.q );

	/* Store our bias and converted angular rates */
	this->accel	= accel;
	this->bias	= x.bias;
	this->pqr	= pqr_raw - x.bias;
}


void
AHRS::compass_update(
	double			heading
)
{
	const state_t &		x( this->state );

	/* Compute the DCM and angle estimate here */
	const Matrix<3,3> 	DCM( quatDC( x.q ) );
	const Vector<3>		THETAe( quat2euler( x.q ) );

	this->kalman_compass_update(
		heading,
		DCM,
		THETAe
	);
}


}
