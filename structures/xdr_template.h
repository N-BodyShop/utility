/** @file xdr_template.h
 Provides a inline templated function xdr_template() which performs
 XDR conversion of a generic value.  Numerous specializations
 of this inline template are provided for many common types, including
 Vector3D and OrientedBox (of arbitrary types themselves).  If
 you attempt to use this function for a type that does not have
 a provided specialization, you will get a link error complaining
 about the function xdr_template() for the type you tried to use.
 You can provide specializations for your own types by copying the
 syntax used below.
 This function allows you to call xdr_template() whenever you wish
 to do XDR conversion, and not worry about the underlying calls
 to xdr_int, xdr_float, etc.
 */

#ifndef XDR_TEMPLATE_H
#define XDR_TEMPLATE_H

#include <rpc/rpc.h>

#include "Vector3D.h"
#include "OrientedBox.h"

template <typename T>
inline bool_t xdr_template(XDR* xdrs, T* val);

template <>
inline bool_t xdr_template(XDR* xdrs, unsigned char* val) {
	return xdr_u_char(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, char* val) {
	return xdr_char(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, unsigned short* val) {
	return xdr_u_short(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, short* val) {
	return xdr_short(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, unsigned int* val) {
	return xdr_u_int(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, int* val) {
	return xdr_int(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, u_int64_t* val) {
	return xdr_u_hyper(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, int64_t* val) {
	return xdr_hyper(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, float* val) {
	return xdr_float(xdrs, val);
}

template <>
inline bool_t xdr_template(XDR* xdrs, double* val) {
	return xdr_double(xdrs, val);
}

template <typename T>
inline bool_t xdr_template(XDR* xdrs, Vector3D<T>* val) {
	return (xdr_template(xdrs, &(val->x))
			&& xdr_template(xdrs, &(val->y))
			&& xdr_template(xdrs, &(val->z)));
}

template <typename T>
inline bool_t xdr_template(XDR* xdrs, OrientedBox<T>* val) {
	return (xdr_template(xdrs, &(val->lesser_corner))
			&& xdr_template(xdrs, &(val->greater_corner)));
}

#endif //XDR_TEMPLATE_H
