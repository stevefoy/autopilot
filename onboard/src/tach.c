/* -*- indent-tabs-mode:T; c-basic-offset:8; tab-width:8; -*- vi: set ts=8:
 *  $Id: tach.c,v 2.0 2002/09/22 02:10:18 tramm Exp $
 *
 * Watches an open-collector Hall effect sensor on PB0 and counts
 * the number of times that it changes per second.
 *
 *************
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

#include "avr.h"
#include "pulse.h"

int main( void )
{
	avr_init();
	pulse_init();

	puts( "$Id: tach.c,v 2.0 2002/09/22 02:10:18 tramm Exp $" );
	putnl();

	while( 1 )
	{
		puthex( high_bits );
		puthex( time() );
		putc( ':' );

		puthexs( pulse_0 );
		putc( ' ' );
		puthexs( pulse_1 );
		putnl();

		pulse_1 = pulse_0 = 0;

		msleep( 8192 );
	}

	return 0;
}
