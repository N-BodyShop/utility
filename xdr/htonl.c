/* Copyright (C) 1993, 1997 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with the GNU C Library; see the file COPYING.LIB.  If not,
   write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

#include <stdint.h>
#include <netinet/in.h>

#undef	htonl
#undef	ntohl

u_int32_t
htonl (x)
     u_int32_t x;
{
#if BYTE_ORDER == BIG_ENDIAN
  return x;
#elif BYTE_ORDER == LITTLE_ENDIAN
  return __bswap_32 (x);
#else
# error "What kind of system is this?"
#endif
}

u_int32_t
ntohl (x)
     u_int32_t x;
{
#if BYTE_ORDER == BIG_ENDIAN
  return x;
#elif BYTE_ORDER == LITTLE_ENDIAN
  return __bswap_32 (x);
#else
# error "What kind of system is this?"
#endif
}
