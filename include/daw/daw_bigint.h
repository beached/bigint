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

#include <daw/daw_algorithm.h>
#include <daw/daw_bit.h>
#include <daw/daw_bounded_string.h>
#include <daw/daw_bounded_vector.h>
#include <daw/daw_exception.h>
#include <daw/daw_math.h>
#include <daw/daw_parser_helper_sv.h>
#include <daw/daw_span.h>
#include <daw/daw_string_view.h>
#include <daw/daw_traits.h>

#include "bigint_impl/bigint_impl.h"

namespace daw::bigint {
	template<size_t BitsNeeded>
	struct bigint_t {
		using value_t = bigint_impl::half_max_t<uintmax_t>;

		static_assert( not std::is_signed_v<value_t>,
		               "Unsupported T, must be unsigned" );

		static_assert( sizeof( value_t ) * 2 <= sizeof( uintmax_t ),
		               "T multiplied by a T must fit into a uintmax_t" );

	private:
		static inline constexpr size_t const m_capacity =
		  bigint_impl::elements_needed<value_t, BitsNeeded>( );

		bigint_impl::bigint_storage_t<value_t, m_capacity + 1> m_data{};

		constexpr void clear_mdata( ) noexcept {
			m_data = decltype( m_data ){};
		}

	public:
		static inline constexpr size_t const capacity_bits = m_capacity * 8;

		[[nodiscard]] constexpr bool is_zero( ) const noexcept {
			return m_data.empty( ) or
			       ( m_data.size( ) = 1 and m_data.front( ) == 0UL );
		}

		constexpr void negate( ) noexcept {
			if( not is_zero( ) ) {
				m_data.sign_flip( );
			}
		}

		constexpr void operator-( ) noexcept {
			return negate( );
		}

		[[nodiscard]] constexpr bigint_t operator~( ) const noexcept {
			bigint_t result{};
			for( auto item : m_data ) {
				result.push_back( ~item );
			}
			for( size_t n = m_data.size( ); n < m_data.capacity( ); ++n ) {
				result.push_back( std::numeric_limits<value_t>::max( ) );
			}
			return result;
		}

		constexpr bigint_t( ) noexcept {
			m_data.push_back( 0 );
		}

		constexpr bool is_negative( ) const noexcept {
			return m_data.m_sign == sign_t::negative;
		}

		template<typename SignedInteger,
		         std::enable_if_t<std::is_signed_v<remove_cvref_t<SignedInteger>>,
		                          std::nullptr_t> = nullptr>
		explicit constexpr bigint_t( SignedInteger v ) {

			m_data.m_sign = v < 0 ? sign_t::negative : sign_t::positive;
			uintmax_t value = 0;
			if( v == std::numeric_limits<intmax_t>::min( ) ) {
				// On intmax_t one cannot multiple min( ) by -1 as max( ) is abs( min(
				// ) ) - 1
				value =
				  static_cast<uintmax_t>( std::numeric_limits<intmax_t>::max( ) ) +
				  1ULL;
			} else {
				value =
				  static_cast<uintmax_t>( v * static_cast<intmax_t>( m_data.m_sign ) );
			}

			size_t elem_needed =
			  bigint_impl::rdiv( bsizeof<SignedInteger>, bsizeof<value_t> );

			while( elem_needed > 0 ) {
				--elem_needed;
				m_data.push_back( bigint_impl::overflow( value ) );
			}
			daw::exception::dbg_precondition_check( value == 0 );
		}

		template<
		  typename UnsignedInteger,
		  std::enable_if_t<
		    all_true_v<std::is_integral_v<remove_cvref_t<UnsignedInteger>>,
		               not std::is_signed_v<remove_cvref_t<UnsignedInteger>>>,
		    std::nullptr_t> = nullptr>
		explicit constexpr bigint_t( UnsignedInteger v ) {

			uintmax_t value = v;
			m_data.m_sign = value < 0 ? sign_t::negative : sign_t::positive;
			value *= static_cast<uintmax_t>( m_data.m_sign );

			while( value > 0 ) {
				m_data.push_back( bigint_impl::overflow( value ) );
			}
			daw::exception::dbg_precondition_check( value == 0 );
		}

		// TODO: add octal/hex/binary via 0o, 0x, 0b
		template<typename CharT>
		explicit constexpr bigint_t( basic_string_view<CharT> str ) {

			str = daw::parser::trim_left( str );
			daw::exception::precondition_check( not str.empty( ) );
			if( str.front( ) == '-' ) {
				str.remove_prefix( );
				m_data.m_sign = sign_t::negative;
			} else {
				if( str.front( ) == '+' ) {
					str.remove_prefix( );
				}
				m_data.m_sign = sign_t::positive;
			}

			// TODO non char input and other bases
			while( not str.empty( ) and daw::parser::is_number( str.front( ) ) ) {
				auto digit =
				  static_cast<uintmax_t>( bigint_impl::to_digit( str.pop_front( ) ) );
				bigint_impl::mul( m_data, 10UL );
				bigint_impl::add( m_data, digit );
			}
		}

		template<typename CharT, size_t N>
		explicit constexpr bigint_t( CharT const ( &str )[N] )
		  : bigint_t(
		      daw::basic_string_view<CharT>( str, str[N - 1] == 0 ? N - 1 : N ) ) {}

		explicit constexpr operator intmax_t( ) const {

			daw::exception::precondition_check(
			  m_data.empty( ) or
			  ( m_data.size( ) - 1U ) * sizeof( value_t ) < sizeof( intmax_t ) );

			intmax_t result = 0;

			auto pos = static_cast<intmax_t>( m_data.size( ) );
			while( --pos >= 0 ) {
				result <<= bsizeof<value_t>;
				result |= m_data[pos];
			}
			// TODO make sure this is not going to be intmax_t::min( )
			if( m_data.m_sign == sign_t::negative ) {
				result = -result;
			}
			return result;
		}

		[[nodiscard]] static constexpr size_t capacity( ) noexcept {
			return m_capacity;
		}

		[[nodiscard]] constexpr size_t size( ) const noexcept {
			return m_data.size( );
		}

		[[nodiscard]] constexpr value_t const &operator[]( size_t idx ) const
		  noexcept {
			return m_data[idx];
		}

		template<size_t>
		friend struct bigint_t;

		template<typename Integer>
		constexpr auto operator*=( Integer &&value )
		  -> std::enable_if_t<std::is_integral_v<remove_cvref_t<Integer>>,
		                      bigint_t &> {

			bigint_t const tmp( std::forward<Integer>( value ) );
			bigint_impl::mul( m_data, tmp.m_data );
			return *this;
		}

		template<typename Integer>
		[[nodiscard]] constexpr auto operator*( Integer &&value ) const
		  -> std::enable_if_t<std::is_integral_v<remove_cvref_t<Integer>>,
		                      bigint_t> {

			auto result = *this;
			auto tmp = bigint_t( std::forward<Integer>( value ) );
			bigint_impl::mul( result.m_data, tmp.m_data );
			return result;
		}

		template<typename Integer>
		constexpr auto operator+=( Integer value )
		  -> std::enable_if_t<std::is_integral_v<Integer>, bigint_t &> {

			bigint_impl::add( m_data,
			                  bigint_t( std::forward<Integer>( value ) ).m_data );
			return *this;
		}

		template<typename Integer>
		[[nodiscard]] constexpr auto operator+( Integer value ) const
		  -> std::enable_if_t<std::is_integral_v<Integer>, bigint_t> {

			auto result = bigint_t( value );
			bigint_impl::add( result.m_data, m_data );
			return result;
		}

		[[nodiscard]] constexpr bigint_t operator+( bigint_t const &value ) const {
			auto result = *this;
			bigint_impl::add( result.m_data, value.m_data );
			return result;
		}

		constexpr bigint_t &set_bit( size_t n ) {
			auto const idx = n / bsizeof<value_t>;
			n -= idx * bsizeof<value_t>;
			for( size_t m = m_data.size( ); m <= idx; ++m ) {
				m_data.push_back( 0 );
			}
			m_data[idx] = daw::set_bits( m_data[idx], n );
			return *this;
		}

		[[nodiscard]] static constexpr bigint_t pow2( size_t n ) {
			bigint_t result{0};
			result.set_bit( n );
			return result;
		}

		[[nodiscard]] static constexpr bigint_t one_shl_minus1( size_t n ) {
			bigint_t result{};
			result.clear_mdata( );
			size_t bit_pos = 0;
			while( bit_pos + bsizeof<value_t> <= n ) {
				result.m_data.push_back( std::numeric_limits<value_t>::max( ) );
				bit_pos += bsizeof<value_t>;
			}
			if( bit_pos == 0 ) {
				result.m_data.push_back( 1 );
				++bit_pos;
			}
			for( ; bit_pos <= n; ++bit_pos ) {
				result.set_bit( bit_pos );
			}
			--result.m_data[0];
			return result;
		}

		template<size_t B>
		[[nodiscard]] constexpr int compare( bigint_t<B> const &rhs ) const
		  noexcept {
			return m_data.compare( rhs.m_data );
		}

		template<size_t B>
		[[nodiscard]] constexpr int compare( bigint_t<B> const &&rhs ) const
		  noexcept {
			return m_data.compare( rhs.m_data );
		}

		template<typename Integer>
		[[nodiscard]] constexpr auto compare( Integer &&rhs ) const noexcept {
			return m_data.compare( bigint_t( std::forward<Integer>( rhs ) ).m_data );
		}

	private:
		[[nodiscard]] static constexpr char nibble_to_char( unsigned n ) {
			if( n < 10 ) {
				return '0' + static_cast<char>( n );
			}
			return 'A' + static_cast<char>( n - 10U );
		}

	public:
		[[nodiscard]] constexpr auto to_hex_string( ) const noexcept {
			using result_t = daw::basic_bounded_string<
			  char,
			  ( static_cast<size_t>( static_cast<double>( BitsNeeded + 4U ) / 4.0 ) +
			    2U )>;
			if( m_data.empty( ) ) {
				return result_t( "0" );
			}
			auto result = result_t( );
			if( is_negative( ) ) {
				result.push_back( '-' );
			}
			bool first_item = true;
			auto const process_nibble = [&]( unsigned nibble ) {
				if( nibble == 0 and first_item ) {
					return;
				}
				first_item = false;
				result.push_back( nibble_to_char( nibble ) );
			};
			for( auto it = m_data.rbegin( ); it != m_data.rend( ); ++it ) {
				auto const cur_v = *it;
				auto cur_q = bigint_impl::int_queue( cur_v );
				if( not first_item ) {
					for( size_t n = cur_q.size( ); n < sizeof( value_t ); ++n ) {
						result.push_back( 0 );
					}
				}
				while( not cur_q.empty( ) ) {
					auto const c = static_cast<unsigned>( cur_q.pop_byte_front( ) );
					// high nibble
					unsigned const n0 = static_cast<unsigned>( c & 0xF0U ) >> 4U;
					process_nibble( n0 );
					// low nibble
					unsigned const n1 = c & 0xFU;
					process_nibble( n1 );
				}
			}
			return result;
		}

		[[nodiscard]] constexpr auto to_decimal_string( ) const noexcept {
			constexpr auto max_digits = base10_digits_in_bits( BitsNeeded );
			using result_t = daw::basic_bounded_string<char, max_digits * 2U>;
			if( m_data.empty( ) ) {
				return result_t( "0" );
			}
			// TODO: When in C++20 use std::vector as stack
			auto int_values = daw::bounded_vector_t<unsigned, max_digits + 1U>( );
			int_values.push_back( 0 );
			bool first_item = true;
			auto const process_nibble = [&]( unsigned nibble ) {
				if( nibble == 0 and first_item ) {
					return;
				}
				first_item = false;
				for( auto &v : int_values ) {
					auto const val = ( v * 16U ) + nibble;
					v = val % 10U;
					nibble = val / 10U;
				}
				while( nibble > 0U ) {
					int_values.push_back( nibble % 10U );
					nibble /= 10U;
				}
			};
			for( auto it = m_data.rbegin( ); it != m_data.rend( ); ++it ) {
				auto const cur_v = *it;
				auto cur_q = bigint_impl::int_queue( cur_v );
				if( not first_item ) {
					for( size_t n = cur_q.size( ); n < sizeof( value_t ); ++n ) {
						process_nibble( 0 );
					}
				}
				while( not cur_q.empty( ) ) {
					auto const c = static_cast<unsigned>( cur_q.pop_byte_front( ) );
					// high nibble
					unsigned const n0 = static_cast<unsigned>( c & 0xF0U ) >> 4U;
					process_nibble( n0 );
					// low nibble
					unsigned const n1 = c & 0xFU;
					process_nibble( n1 );
				}
			}
			auto result = result_t( );
			if( is_negative( ) ) {
				result.push_back( '-' );
			}

			daw::algorithm::transform(
			  int_values.rbegin( ), int_values.rend( ), daw::back_inserter( result ),
			  []( unsigned u ) -> char { return '0' + static_cast<char>( u ); } );
			return result;
		}
	}; // namespace daw::bigint

	template<typename CharT, size_t N>
	bigint_t( CharT const ( & )[N] )
	  ->bigint_t<bigint_impl::bits_needed<N - 1>( )>;

	template<size_t base10_digits>
	using bigint_digits_t = bigint_t<bits_needed_for_digits( base10_digits )>;

	namespace bigint_impl {
		template<size_t N>
		[[nodiscard]] constexpr std::true_type
		is_bigint_test( bigint_t<N> const & ) noexcept;

		template<size_t N>
		[[nodiscard]] constexpr std::true_type
		is_bigint_test( bigint_t<N> && ) noexcept;

		[[nodiscard]] constexpr std::false_type is_bigint_test( ... ) noexcept;
	} // namespace bigint_impl

	template<typename T>
	using is_bigint_t = decltype(
	  bigint_impl::is_bigint_test( std::declval<remove_cvref_t<T>>( ) ) );

	template<typename T>
	inline constexpr bool const is_bigint_v = is_bigint_t<T>::value;

	template<size_t LhsB, size_t RhsB>
	[[nodiscard]] constexpr auto
	operator==( bigint_t<LhsB> const &lhs, bigint_t<RhsB> const &rhs ) noexcept {
		return lhs.compare( rhs ) == 0;
	}

	template<size_t LhsB, typename Integer>
	[[nodiscard]] constexpr auto operator==( bigint_t<LhsB> const &lhs,
	                                         Integer &&rhs ) noexcept
	  -> std::enable_if_t<std::is_integral_v<Integer>, bool> {

		return lhs.compare( std::forward<Integer>( rhs ) ) == 0;
	}

	template<typename Integer, size_t RhsB>
	[[nodiscard]] constexpr auto operator==( Integer &&lhs,
	                                         bigint_t<RhsB> const &rhs ) noexcept
	  -> std::enable_if_t<std::is_integral_v<Integer>, bool> {

		return bigint_t<bsizeof<Integer>>( std::forward<Integer>( lhs ) )
		         .compare( rhs ) == 0;
	}

	template<size_t LhsB, size_t RhsB>
	[[nodiscard]] constexpr auto
	operator!=( bigint_t<LhsB> const &lhs, bigint_t<RhsB> const &rhs ) noexcept {

		return lhs.m_data.compare( rhs ) != 0;
	}

	template<size_t LhsB, typename Integer>
	[[nodiscard]] constexpr auto operator!=( bigint_t<LhsB> const &lhs,
	                                         Integer &&rhs ) noexcept
	  -> std::enable_if_t<std::is_integral_v<Integer>, bool> {
		return lhs.compare( std::forward<Integer>( rhs ) ) != 0;
	}

	template<typename Integer, size_t RhsB>
	[[nodiscard]] constexpr auto operator!=( Integer &&lhs,
	                                         bigint_t<RhsB> const &rhs ) noexcept
	  -> std::enable_if_t<std::is_integral_v<Integer>, bool> {
		return std::forward<Integer>( lhs ).compare(
		         std::forward<Integer>( rhs ) ) != 0;
	}
} // namespace daw::bigint
