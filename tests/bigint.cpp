// The MIT License (MIT)
//
// Copyright (c) 2019 Darrell Wright
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files( the "Software" ), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and / or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <iostream>

#include "daw/daw_bigint.h"

template<size_t Bits>
static void display( daw::bigint::bigint_t<Bits> const &v ) {
	std::cout << "bits: " << Bits;
	auto const ds = v.to_decimal_string( );
	std::cout << " | decimal: " << ds;
	auto const hs = v.to_hex_string( );
	std::cout << " | hex: " << hs;
	std::cout << '\n';
}

int main( ) {
	using daw::bigint::bigint_t;
	constexpr auto t1 = bigint_t( "12345678459485983745983475398574897" );
	display( t1 );
	constexpr auto t2 = t1 + 5678ULL;
	display( t2 );
	constexpr auto t3 = t2;
	display( t3 );
	constexpr auto t4 = t3 + t3;
	display( t4 );
	display( bigint_t( "123456789" ) );
	display( bigint_t( "12345678901" ) );
	display( bigint_t( "123456789012" ) );
	display( bigint_t( "1234567890123" ) );
	display( bigint_t( "12345678901234" ) );
	display( bigint_t( "123456789012345" ) );
	display( bigint_t( "-4294967296" ) );
	display( bigint_t( "1" ) );
	display( bigint_t( "-1" ) );
	display( bigint_t( "-128" ) );
	display(
	  bigint_t( "1234567890123451111111111111111111111111111111111111111111111111"
	            "1111111111111111111111111111000000000000000000000000000000000000"
	            "0000000000000000000000000000000000000000000000000000000" ) );
}
