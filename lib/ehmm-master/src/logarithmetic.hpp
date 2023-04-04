/* Copyright (C) 2010 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef _LOGARITHMETIC_H_
#define _LOGARITHMETIC_H_

#include <cstdlib>
#include <cmath>
#include <limits>

/* Log Sum of Exponentials Algorithm */

template <class RealType = double>
RealType logadd(RealType a, RealType b)
{
        if (a < b) return a == -std::numeric_limits<RealType>::infinity() ? b : b + std::log1p(std::exp(a-b));
        else       return b == -std::numeric_limits<RealType>::infinity() ? a : a + std::log1p(std::exp(b-a));
}

template <class RealType = double>
RealType logsub(RealType a, RealType b)
{
        return b == -std::numeric_limits<RealType>::infinity() ? a : a + std::log(1-exp(b-a));
}

#endif /* _LOGARITHMETIC_H_ */
