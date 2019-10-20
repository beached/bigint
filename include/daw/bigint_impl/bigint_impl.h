// The MIT License (MIT)
//
// Copyright (c) 2018-2019 Darrell Wright
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

#pragma once

#include <array>
#include <iterator>

#include <daw/daw_algorithm.h>
#include <daw/daw_bit.h>
#include <daw/daw_exception.h>
#include <daw/daw_math.h>
#include <daw/daw_move.h>
#include <daw/daw_parser_helper_sv.h>
#include <daw/daw_span.h>
#include <daw/daw_string_view.h>
#include <daw/daw_traits.h>

namespace daw::bigint {
	enum class sign_t : int_fast8_t { positive = 1, negative = -1 };

	// It is used to calculate how many bits are needed to represent a
	// base 10 number with value digits
	constexpr size_t bits_needed_for_digits( size_t value ) noexcept {
		constexpr double const log2_of_10 = 3.32192809488736218;
		// result = log2( 10^value ) + 1;
		return static_cast<size_t>( static_cast<double>( value ) * log2_of_10 ) +
		       1U;
	}

	constexpr size_t base10_digits_in_bits( size_t num_bits ) noexcept {
		constexpr double const log10_of_2 = static_cast<double>(
		  0.301029995663981195213738894724493026768189881462108541310 );
		return static_cast<size_t>( static_cast<double>( num_bits ) * log10_of_2 ) +
		       1U;
	}
} // namespace daw::bigint

namespace daw::bigint::bigint_impl {

	template<typename T>
	class int_queue {
		static_assert( std::is_integral_v<T>,
		               "T must be an unsigned integer type" );
		static_assert( std::is_unsigned_v<T>,
		               "T must be an unsigned integer type" );
		static_assert( sizeof( T ) >= sizeof( unsigned ),
		               "T must be at least as large as an unsigned" );

		T m_value = 0;
		unsigned m_size = 0;

		[[nodiscard]] static constexpr T make_mask( unsigned pos ) noexcept {
			return 0xFFU << ( pos * 8U );
		}
		[[nodiscard]] static constexpr unsigned calc_size( T value ) noexcept {
			if( value == 0 ) {
				return sizeof( T );
			}
			unsigned result = sizeof( T ) - 1U;
			while( result > 0 ) {
				T const mask = make_mask( result );
				T const remain = mask & value;
				if( remain != 0 ) {
					// TODO: see if we have to use CHAR_BITS
					return result + 1U;
				}
				--result;
			}
			return result + 1U;
		}

	public:
		constexpr int_queue( ) noexcept = default;
		constexpr int_queue( T value ) noexcept
		  : m_value( value )
		  , m_size( calc_size( value ) ) {}

		[[nodiscard]] constexpr uint8_t pop_byte_front( ) noexcept {
			T const mask = make_mask( m_size - 1U );
			T result = m_value & mask;
			m_value &= ~mask;
			--m_size;
			result >>= m_size * 8U;
			return static_cast<uint8_t>( result );
		}

		[[nodiscard]] constexpr bool empty( ) const noexcept {
			return m_size < 1U;
		}

		[[nodiscard]] constexpr size_t size( ) const noexcept {
			return m_size;
		}
	};

	static_assert( sizeof( uintmax_t ) == sizeof( uint64_t ),
	               "Assumed that uint64_t is largest unsigned integer type" );
	static_assert( sizeof( intmax_t ) == sizeof( int64_t ),
	               "Assumed that uint64_t is largest integer type" );

	template<typename T = uint8_t>
	[[nodiscard]] constexpr T to_digit( char c ) noexcept {
		if( '0' <= c && c <= '9' ) {
			return static_cast<T>( c - '0' );
		}
		if( 'a' <= c && c <= 'f' ) {
			return static_cast<T>( 10 ) + static_cast<T>( c - 'a' );
		}
		if( 'A' <= c && c <= 'F' ) {
			return static_cast<T>( 10 ) + static_cast<T>( c - 'A' );
		} else { // Work around constexpr bug in gcc 8.2.  Need to put else clause
			std::terminate( );
		}
	}

	template<typename T = uint8_t>
	[[nodiscard]] constexpr T to_digit( wchar_t c ) noexcept {
		if( L'0' <= c && c <= L'9' ) {
			return static_cast<T>( c - L'0' );
		}
		if( L'a' <= c && c <= L'f' ) {
			return static_cast<T>( 10 ) + static_cast<T>( c - L'a' );
		}
		if( 'A' <= c && c <= 'F' ) {
			return static_cast<T>( 10 ) + static_cast<T>( c - L'A' );
		} else { // Work around constexpr bug in gcc 8.2.  Need to put else clause
			std::terminate( );
		}
	}

	template<typename T>
	using half_max_t = std::conditional_t<
	  sizeof( T ) == 8, uint32_t,
	  std::conditional_t<sizeof( T ) == 4, uint16_t,
	                     std::conditional_t<sizeof( T ) == 2, uint8_t, void>>>;

	[[nodiscard]] constexpr auto low_part( uintmax_t value ) noexcept {
		using value_t = half_max_t<uintmax_t>;
		auto const lower_mask =
		  daw::mask_msb<uintmax_t>( bsizeof<uintmax_t> - bsizeof<value_t> );

		return static_cast<value_t>( value bitand lower_mask );
	}

	[[nodiscard]] constexpr auto high_part( uintmax_t value ) noexcept {
		using value_t = half_max_t<uintmax_t>;
		auto const high_mask =
		  daw::mask_lsb<uintmax_t>( bsizeof<uintmax_t> - bsizeof<value_t> );

		return static_cast<value_t>( value bitand high_mask );
	}

	template<typename T, size_t ItemCount>
	struct bigint_storage_t {
		using value_t = T;
		using size_type = size_t;
		using reference = value_t &;
		using const_reference = value_t const &;
		using iterator = typename std::array<value_t, ItemCount>::iterator;
		using const_iterator =
		  typename std::array<value_t, ItemCount>::const_iterator;
		using reverse_iterator = std::reverse_iterator<iterator>;
		using const_reverse_iterator = std::reverse_iterator<const_iterator>;

		static_assert( ItemCount > 0, "Cannot store 0 of anything" );

		std::array<value_t, ItemCount> m_data = {0};
		size_type m_idx = 0;
		sign_t m_sign = sign_t::positive;

		constexpr bigint_storage_t( ) noexcept = default;

		template<size_t N,
		         std::enable_if_t<( N != ItemCount ), std::nullptr_t> = nullptr>
		explicit constexpr bigint_storage_t( bigint_storage_t<T, N> const &other )
		  : m_idx( other.m_idx )
		  , m_sign( other.m_sign ) {
			daw::exception::precondition_check<std::out_of_range>( m_idx <=
			                                                       ItemCount );
			daw::algorithm::copy_n( other.m_data.data( ), m_data.data( ), m_idx );
		}

		template<size_t N,
		         std::enable_if_t<( N != ItemCount ), std::nullptr_t> = nullptr>
		constexpr bigint_storage_t &operator=( bigint_storage_t<T, N> const &rhs ) {
			daw::exception::precondition_check<std::out_of_range>( rhs.m_idx <=
			                                                       ItemCount );
			if( this == &rhs ) {
				return *this;
			}

			m_idx = rhs.m_idx;
			m_sign = rhs.m_sign;
			daw::algorithm::copy_n( rhs.m_data.data( ), m_data.data( ), m_idx );
			return *this;
		}

		template<size_t N,
		         std::enable_if_t<( N != ItemCount ), std::nullptr_t> = nullptr>
		explicit constexpr bigint_storage_t( bigint_storage_t<T, N> &&other )
		  : m_idx( daw::move( other.m_idx ) )
		  , m_sign( daw::move( other.m_sign ) ) {
			daw::exception::precondition_check<std::out_of_range>( m_idx <=
			                                                       ItemCount );

			daw::algorithm::move_n( other.m_data.data( ), m_data.data( ), m_idx );
		}

		template<size_t N,
		         std::enable_if_t<( N != ItemCount ), std::nullptr_t> = nullptr>
		constexpr bigint_storage_t &operator=( bigint_storage_t<T, N> &&rhs ) {
			daw::exception::precondition_check<std::out_of_range>( rhs.m_idx <=
			                                                       ItemCount );
			if( this == &rhs ) {
				return *this;
			}

			m_idx = daw::move( rhs.m_idx );
			m_sign = daw::move( rhs.m_sign );
			daw::algorithm::move_n( rhs.m_data.data( ), m_data.data( ), m_idx );
			return *this;
		}

		constexpr void clear( ) noexcept {
			while( m_idx > 0 ) {
				--m_idx;
				m_data[m_idx] = 0;
			}
		}

		template<typename Size>
		[[nodiscard]] constexpr reference operator[]( Size idx ) noexcept {
			daw::exception::dbg_precondition_check<std::out_of_range>(
			  in_range( static_cast<size_t>( idx ), 0ULL, ItemCount ) );
			return m_data[static_cast<size_type>( idx )];
		}

		template<typename Size>
		[[nodiscard]] constexpr const_reference operator[]( Size idx ) const
		  noexcept {
			daw::exception::dbg_precondition_check<std::out_of_range>(
			  in_range( static_cast<size_t>( idx ), 0ULL, ItemCount ) );
			return m_data[static_cast<size_type>( idx )];
		}

		[[nodiscard]] constexpr size_type const &size( ) const noexcept {
			return m_idx;
		}

		[[nodiscard]] constexpr size_type &size( ) noexcept {
			return m_idx;
		}

		[[nodiscard]] constexpr bool empty( ) const noexcept {
			return m_idx == 0;
		}

		[[nodiscard]] constexpr iterator begin( ) noexcept {
			return m_data.begin( );
		}

		[[nodiscard]] constexpr const_iterator begin( ) const noexcept {
			return m_data.begin( );
		}

		[[nodiscard]] constexpr reverse_iterator rbegin( ) noexcept {
			return reverse_iterator( end( ) );
		}

		[[nodiscard]] constexpr const_reverse_iterator rbegin( ) const noexcept {
			return const_reverse_iterator( end( ) );
		}

		[[nodiscard]] constexpr iterator end( ) noexcept {
			return m_data.begin( ) + static_cast<ptrdiff_t>( m_idx );
		}

		[[nodiscard]] constexpr const_iterator end( ) const noexcept {
			return m_data.begin( ) + static_cast<ptrdiff_t>( m_idx );
		}

		[[nodiscard]] constexpr reverse_iterator rend( ) noexcept {
			return reverse_iterator( begin( ) );
		}

		[[nodiscard]] constexpr const_reverse_iterator rend( ) const noexcept {
			return const_reverse_iterator( begin( ) );
		}

		[[nodiscard]] constexpr size_type capacity( ) const noexcept {
			return ItemCount;
		}

		[[nodiscard]] constexpr bool full( ) noexcept {
			return size( ) >= capacity( );
		}

		constexpr void sign_flip( ) noexcept {
			m_sign = m_sign == sign_t::positive ? sign_t::negative : sign_t::positive;
		}

		template<typename Value, size_t Brhs>
		[[nodiscard]] constexpr int
		unsigned_compare( bigint_storage_t<Value, Brhs> const &rhs ) const {
			if( size( ) > rhs.size( ) ) {
				return 1;
			} else if( size( ) < rhs.size( ) ) {
				return -1;
			}
			for( size_t n = size( ); n > 0; --n ) {
				if( m_data[n] > rhs.m_data[n] ) {
					return 1;
				} else if( m_data[n] < rhs.m_data[n] ) {
					return -1;
				}
			}
			return 0;
		}

		template<typename Value, size_t Brhs>
		[[nodiscard]] constexpr int
		compare( bigint_storage_t<Value, Brhs> const &rhs ) const {
			if( m_sign == sign_t::positive ) {
				if( rhs.m_sign == sign_t::negative ) {
					// pos neg
					return 1;
				}
				// pos pos
				return unsigned_compare( rhs );
			} else if( rhs.m_sign == sign_t::positive ) {
				// neg pos
				return -1;
			}
			// neg neg
			return -unsigned_compare( rhs );
		}

		constexpr void push_back( value_t v ) {
			daw::exception::dbg_precondition_check<std::overflow_error>(
			  not full( ) );

			m_data[m_idx++] = v;
		}

		[[nodiscard]] constexpr reference back( ) {
			daw::exception::dbg_precondition_check<std::out_of_range>(
			  !empty( ), "Attempt to access value on empty stack" );
			return m_data[m_idx - 1];
		}

		[[nodiscard]] constexpr const_reference back( ) const {
			daw::exception::dbg_precondition_check<std::out_of_range>(
			  !empty( ), "Attempt to access value on empty stack" );
			return m_data[m_idx - 1];
		}

		[[nodiscard]] constexpr value_t pop_back( ) {
			daw::exception::dbg_precondition_check<std::overflow_error>(
			  !empty( ), "Attempt to pop empty stack" );
			auto result = back( );
			back( ) = 0;
			--m_idx;
			return result;
		}
	};

	template<typename T, size_t BitsNeeded>
	[[nodiscard]] constexpr size_t elements_needed( ) noexcept {
		size_t const bits_in_T = bsizeof<T>;
		size_t result = ( BitsNeeded + ( bits_in_T / 2 ) ) / bits_in_T;
		if( result == 0 ) {
			return 1;
		}
		return result;
	}

	template<size_t digits>
	[[nodiscard]] constexpr size_t bits_needed( ) noexcept {
		using value_t = half_max_t<uintmax_t>;
		constexpr size_t num_bits = daw::bigint::bits_needed_for_digits( digits );
		return elements_needed<value_t, num_bits>( ) * sizeof( value_t ) * 8U;
	}

	template<typename T, typename U>
	[[nodiscard]] constexpr auto rdiv( T lhs, U rhs ) noexcept {
		auto result = lhs / rhs;
		if( result * rhs < lhs ) {
			++result;
		}
		return result;
	}

	[[nodiscard]] constexpr auto overflow( uintmax_t &carry ) noexcept {
		using value_t = half_max_t<uintmax_t>;
		auto result = low_part( carry );
		carry >>= bsizeof<value_t>;
		return result;
	}

	template<typename value_t, size_t N>
	constexpr void add( bigint_storage_t<value_t, N> &lhs, value_t rhs,
	                    size_t index = 0 ) {

		daw::exception::precondition_check( rhs >= 0 );

		auto carry = static_cast<uintmax_t>( rhs );

		for( ; index < lhs.size( ); ++index ) {
			carry += lhs[index];
			lhs[index] = overflow( carry );
		}
		if( carry > 0 ) {
			lhs.push_back( static_cast<value_t>( carry ) );
		}
	}

	template<typename value_t, size_t N>
	[[nodiscard]] constexpr uintmax_t div( bigint_storage_t<value_t, N> &lhs,
	                                       value_t const rhs ) {
		if( lhs.m_data.size( ) == 0 ) {
			return rhs;
		}
		bigint_storage_t<value_t, N> result{};
		uint64_t carry = 0;
		for( size_t n = lhs.m_data.size( ); n > 0; --n ) {
			carry += lhs.m_data[n - 1];
			auto tmp = carry / rhs;
			result.push_back( tmp );
			carry -= tmp * rhs;
		}
		lhs.m_data = daw::move( result );
		return carry;
	};

	template<typename value_t, size_t LhsN, size_t RhsN>
	[[nodiscard]] constexpr bigint_storage_t<value_t, LhsN>
	div( bigint_storage_t<value_t, LhsN> &lhs,
	     bigint_storage_t<value_t, RhsN> const &rhs ) {

		// result = lhs/rhs
		// lhs /= rhs;
		auto const cmp = lhs.compare( rhs );
		if( cmp < 0 ) {
			lhs = bigint_storage_t<value_t, LhsN>{};
			return rhs;
		} else if( cmp == 0 ) {
			lhs = 1UL;
			return 0UL;
		}
		// lhs > rhs
		daw::exception::precondition_check( rhs != 0 );
		daw::span<value_t> dividend( lhs.front( ), lhs.size( ) );
		daw::span<value_t> divisor( rhs.front( ), rhs.size( ) );

		while( !dividend.empty( ) ) {
			auto cur_quotient = dividend.back( ) / divisor.back( );
		}
	};

	template<typename value_t, size_t N>
	constexpr void sub( bigint_storage_t<value_t, N> &lhs, value_t ) {
		Unused( lhs ); // TODO Fix this mess
		               /*
		                     auto const lower_mask =
		                         daw::get_left_mask<uintmax_t>( bsizeof<uintmax_t> -
		                  bsizeof<value_t> );
		           
		                     if( lhs_idx == 0 ) {
		                       daw::exception::precondition_check( rhs == 0, "We should never be
		                  calling sub with a rhs > lhs and lhs == 0" ); return;
		                     }
		                     if( static_cast<size_t>( bit_pos ) <= lhs_idx ) {
		           
		                     }
		                     uintmax_t carry = rhs;
		                     for( size_t n=static_cast<size_t>( bit_pos ); n<idx; ++n ) {
		                       if( lhs[n] >= carry ) {
		                         lhs[n] -= carry;
		                       } else {
		           
		                       }
		                       carry += lhs[n];
		                       lhs[n] = static_cast<value_t>( carry bitand lower_mask );
		                       carry >>= bsizeof<value_t>;
		                     }
		                     if( carry > 0 ) {
		                       daw::exception::precondition_check<std::overflow_error>( idx <
		                  lhs.size( ) ); lhs[idx++] = static_cast<value_t>( carry );
		                     }
		                     */
	}

	template<typename value_t, size_t N>
	constexpr void add( bigint_storage_t<value_t, N> &lhs, uintmax_t rhs ) {

		auto const lower = low_part( rhs );
		if( lower > 0 ) {
			add( lhs, lower, 0 );
		}

		auto const higher = high_part( rhs );
		if( higher > 0 ) {
			add( lhs, higher, 1 );
		}
	}

	template<typename value_t, size_t N>
	constexpr void add( bigint_storage_t<value_t, N> &lhs,
	                    bigint_storage_t<value_t, N> const &rhs ) {

		for( size_t n = 0; n < rhs.size( ); ++n ) {
			add( lhs, rhs[n], n );
		}
	}

	template<typename value_t, size_t N>
	constexpr void mul( bigint_storage_t<value_t, N> &lhs, value_t rhs,
	                    size_t index ) {

		if( rhs == 0ULL ) {
			lhs.clear( );
			return;
		}

		uintmax_t carry = 0;

		auto tmp = lhs;
		for( size_t pos = 0ULL; pos < lhs.size( ); ++pos ) {
			carry +=
			  static_cast<uintmax_t>( lhs[pos] ) * static_cast<uintmax_t>( rhs );
			tmp[pos] = overflow( carry );
		}
		if( carry > 0 ) {
			tmp.push_back( static_cast<value_t>( carry ) );
			carry = 0;
		}
		lhs = tmp;
		if( index > 0 ) {
			while( index > 0 ) {
				--index;
				for( size_t pos = 0; pos < lhs.size( ); ++pos ) {
					carry +=
					  static_cast<uintmax_t>( lhs[pos] ) *
					  static_cast<uintmax_t>( std::numeric_limits<value_t>::max( ) );
					tmp[pos] = overflow( carry );
				}
			}
			if( carry > 0 ) {
				tmp.push_back( static_cast<value_t>( carry ) );
			}
			lhs = tmp;
		}
	}

	template<typename value_t, size_t N>
	constexpr void mul( bigint_storage_t<value_t, N> &lhs, uintmax_t rhs ) {

		auto tmp = lhs;

		mul( lhs, low_part( rhs ), 0 );
		mul( tmp, high_part( rhs ), 1 );
		add( lhs, tmp );
	}

	template<typename value_t, size_t N>
	constexpr void mul( bigint_storage_t<value_t, N> &lhs,
	                    bigint_storage_t<value_t, N> const &rhs ) {

		if( lhs.empty( ) ) {
			return;
		}
		if( rhs.empty( ) ) {
			lhs.clear( );
			return;
		}

		bigint_storage_t<value_t, N> result{};
		result.push_back( 0 );

		for( size_t n = 0; n < rhs.size( ); ++n ) {
			auto tmp = lhs;
			mul( tmp, rhs[n], n );
			add( result, tmp );
		}
		lhs = daw::move( result );
	}
} // namespace daw::bigint::bigint_impl
