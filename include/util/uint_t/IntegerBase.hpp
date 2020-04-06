/*******************************************************************************
 * include/uint_t/IntegerBase.hpp
 *
 * Copyright (C) 2016 Marvin Löbel
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

namespace pwm {

static_assert(sizeof(int) * 8 == 32, "Make sure the logic here remains correct");

template<class T, class X = void>
struct ConstIntegerBaseTrait {
};

template<class T, class X = void>
struct IntegerBaseTrait: public ConstIntegerBaseTrait<T, X> {
};

template<class Self>
class IntegerBaseWithSelf;
template<class Self, class Other>
class IntegerBaseWith32;
template<class Self, class Other>
class IntegerBaseWith64;

template<class Self, class Other>
class ConstIntegerBaseWith32 {
public:
    typedef typename ConstIntegerBaseTrait<Self>::Dispatch::SelfMaxBit SelfMaxBit;
private:
    inline static SelfMaxBit cast_for_32_op(const Self& self) {
        return ConstIntegerBaseTrait<Self>::Dispatch::template cast_for_op<Self, SelfMaxBit>(self);
    }

    friend class IntegerBaseWith32<Self, Other>;
public:
    friend SelfMaxBit operator+(const Other& lhs, const Self& rhs) { return lhs                   + cast_for_32_op(rhs);   }
    friend SelfMaxBit operator+(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   + rhs;                   }

    friend SelfMaxBit operator-(const Other& lhs, const Self& rhs) { return lhs                   - cast_for_32_op(rhs);   }
    friend SelfMaxBit operator-(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   - rhs;                   }

    friend SelfMaxBit operator*(const Other& lhs, const Self& rhs) { return lhs                   * cast_for_32_op(rhs);   }
    friend SelfMaxBit operator*(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   * rhs;                   }

    friend SelfMaxBit operator/(const Other& lhs, const Self& rhs) { return lhs                   / cast_for_32_op(rhs);   }
    friend SelfMaxBit operator/(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   / rhs;                   }

    friend SelfMaxBit operator%(const Other& lhs, const Self& rhs) { return lhs                   % cast_for_32_op(rhs);   }
    friend SelfMaxBit operator%(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   % rhs;                   }

    friend SelfMaxBit operator&(const Other& lhs, const Self& rhs) { return lhs                   & cast_for_32_op(rhs);   }
    friend SelfMaxBit operator&(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   & rhs;                   }

    friend SelfMaxBit operator|(const Other& lhs, const Self& rhs) { return lhs                   | cast_for_32_op(rhs);   }
    friend SelfMaxBit operator|(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   | rhs;                   }

    friend SelfMaxBit operator^(const Other& lhs, const Self& rhs) { return lhs                   ^ cast_for_32_op(rhs);   }
    friend SelfMaxBit operator^(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   ^ rhs;                   }

    friend SelfMaxBit operator<<(const Other& lhs, const Self& rhs) { return lhs                   << cast_for_32_op(rhs);   }
    friend SelfMaxBit operator<<(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   << rhs;                   }

    friend SelfMaxBit operator>>(const Other& lhs, const Self& rhs) { return lhs                   >> cast_for_32_op(rhs);   }
    friend SelfMaxBit operator>>(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   >> rhs;                   }

    friend bool operator==(const Other& lhs, const Self& rhs) { return lhs                   == cast_for_32_op(rhs);   }
    friend bool operator==(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   == rhs;                   }

    friend bool operator!=(const Other& lhs, const Self& rhs) { return lhs                   != cast_for_32_op(rhs);   }
    friend bool operator!=(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   != rhs;                   }

    friend bool operator>(const Other& lhs, const Self& rhs) { return lhs                   > cast_for_32_op(rhs);   }
    friend bool operator>(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   > rhs;                   }

    friend bool operator<(const Other& lhs, const Self& rhs) { return lhs                   < cast_for_32_op(rhs);   }
    friend bool operator<(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   < rhs;                   }

    friend bool operator>=(const Other& lhs, const Self& rhs) { return lhs                   >= cast_for_32_op(rhs);   }
    friend bool operator>=(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   >= rhs;                   }

    friend bool operator<=(const Other& lhs, const Self& rhs) { return lhs                   <= cast_for_32_op(rhs);   }
    friend bool operator<=(const Self& lhs, const Other& rhs) { return cast_for_32_op(lhs)   <= rhs;                   }
};

template<class Self, class Other>
class ConstIntegerBaseWith64 {
public:
    typedef typename ConstIntegerBaseTrait<Self>::Dispatch::SelfMaxBit SelfMaxBit;
private:
    inline static Other cast_for_64_op(const Self& self) {
        return ConstIntegerBaseTrait<Self>::Dispatch::template cast_for_op<Self, uint64_t>(self);
    }

    friend class IntegerBaseWith64<Self, Other>;
public:
    friend Other   operator+(const Other& lhs, const Self& rhs) { return lhs                   + cast_for_64_op(rhs);   }
    friend Other   operator+(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   + rhs;                   }

    friend Other   operator-(const Other& lhs, const Self& rhs) { return lhs                   - cast_for_64_op(rhs);   }
    friend Other   operator-(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   - rhs;                   }

    friend Other   operator*(const Other& lhs, const Self& rhs) { return lhs                   * cast_for_64_op(rhs);   }
    friend Other   operator*(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   * rhs;                   }

    friend Other   operator/(const Other& lhs, const Self& rhs) { return lhs                   / cast_for_64_op(rhs);   }
    friend Other   operator/(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   / rhs;                   }

    friend Other   operator%(const Other& lhs, const Self& rhs) { return lhs                   % cast_for_64_op(rhs);   }
    friend Other   operator%(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   % rhs;                   }

    friend Other   operator&(const Other& lhs, const Self& rhs) { return lhs                   & cast_for_64_op(rhs);   }
    friend Other   operator&(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   & rhs;                   }

    friend Other   operator|(const Other& lhs, const Self& rhs) { return lhs                   | cast_for_64_op(rhs);   }
    friend Other   operator|(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   | rhs;                   }

    friend Other   operator^(const Other& lhs, const Self& rhs) { return lhs                   ^ cast_for_64_op(rhs);   }
    friend Other   operator^(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   ^ rhs;                   }

    friend Other   operator<<(const Other& lhs, const Self& rhs) { return lhs                   << cast_for_64_op(rhs);   }
    friend Other   operator<<(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   << rhs;                   }

    friend Other   operator>>(const Other& lhs, const Self& rhs) { return lhs                   >> cast_for_64_op(rhs);   }
    friend Other   operator>>(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   >> rhs;                   }

    friend bool operator==(const Other& lhs, const Self& rhs) { return lhs                   == cast_for_64_op(rhs);   }
    friend bool operator==(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   == rhs;                   }

    friend bool operator!=(const Other& lhs, const Self& rhs) { return lhs                   != cast_for_64_op(rhs);   }
    friend bool operator!=(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   != rhs;                   }

    friend bool operator>(const Other& lhs, const Self& rhs) { return lhs                   > cast_for_64_op(rhs);   }
    friend bool operator>(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   > rhs;                   }

    friend bool operator<(const Other& lhs, const Self& rhs) { return lhs                   < cast_for_64_op(rhs);   }
    friend bool operator<(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   < rhs;                   }

    friend bool operator>=(const Other& lhs, const Self& rhs) { return lhs                   >= cast_for_64_op(rhs);   }
    friend bool operator>=(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   >= rhs;                   }

    friend bool operator<=(const Other& lhs, const Self& rhs) { return lhs                   <= cast_for_64_op(rhs);   }
    friend bool operator<=(const Self& lhs, const Other& rhs) { return cast_for_64_op(lhs)   <= rhs;                   }
};

template<class Self>
class ConstIntegerBaseWithSelf {
public:
    typedef typename ConstIntegerBaseTrait<Self>::Dispatch::SelfMaxBit SelfMaxBit;
private:
    inline static SelfMaxBit cast_for_self_op(const Self& self) {
        return ConstIntegerBaseTrait<Self>::Dispatch::template cast_for_op<Self, SelfMaxBit>(self);
    }

    friend class IntegerBaseWithSelf<Self>;
public:
    friend SelfMaxBit operator+(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) + cast_for_self_op(rhs); }
    friend SelfMaxBit operator-(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) - cast_for_self_op(rhs); }
    friend SelfMaxBit operator*(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) * cast_for_self_op(rhs); }
    friend SelfMaxBit operator/(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) / cast_for_self_op(rhs); }
    friend SelfMaxBit operator%(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) % cast_for_self_op(rhs); }
    friend SelfMaxBit operator&(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) & cast_for_self_op(rhs); }
    friend SelfMaxBit operator|(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) | cast_for_self_op(rhs); }
    friend SelfMaxBit operator^(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) ^ cast_for_self_op(rhs); }
    friend SelfMaxBit operator~(const Self& self) { return ~(cast_for_self_op(self)); }
    friend SelfMaxBit operator<<(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) << cast_for_self_op(rhs); }
    friend SelfMaxBit operator>>(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) >> cast_for_self_op(rhs); }
    friend bool operator==(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) == cast_for_self_op(rhs); }
    friend bool operator!=(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) != cast_for_self_op(rhs); }
    friend bool operator>(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) > cast_for_self_op(rhs); }
    friend bool operator<(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) < cast_for_self_op(rhs); }
    friend bool operator>=(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) >= cast_for_self_op(rhs); }
    friend bool operator<=(const Self& lhs, const Self& rhs) { return cast_for_self_op(lhs) <= cast_for_self_op(rhs); }
};

template<typename...>
class ConstIntegerBaseCombiner;

template<typename T, typename... Ts>
class ConstIntegerBaseCombiner<T, Ts...>: public T, public ConstIntegerBaseCombiner<Ts...> {
};

template<typename T>
class ConstIntegerBaseCombiner<T>: public T {
};

template<class Self>
using ConstIntegerBase = ConstIntegerBaseCombiner<
    ConstIntegerBaseWithSelf<Self>,
    ConstIntegerBaseWith32<Self, unsigned char>,
    ConstIntegerBaseWith32<Self, char>,
    ConstIntegerBaseWith32<Self, signed char>,
    ConstIntegerBaseWith32<Self, unsigned short int>,
    ConstIntegerBaseWith32<Self, signed short int>,
    ConstIntegerBaseWith32<Self, unsigned int>,
    ConstIntegerBaseWith32<Self, signed int>,
    ConstIntegerBaseWith64<Self, unsigned long int>,
    ConstIntegerBaseWith64<Self, signed long int>,
    ConstIntegerBaseWith64<Self, unsigned long long int>,
    ConstIntegerBaseWith64<Self, signed long long int>
>;

template<class Self, class Other>
class IntegerBaseWith32: public ConstIntegerBaseWith32<Self, Other> {
private:
    inline static void assign(Self& self, uint32_t v) {
        IntegerBaseTrait<Self>::Dispatch::template assign<Self, uint32_t>(self, v);
    }

    inline static void assign(Self& self, uint64_t v) {
        IntegerBaseTrait<Self>::Dispatch::template assign<Self, uint64_t>(self, v);
    }
public:
    Self& operator+=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self + v); return self; }
    Self& operator-=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self - v); return self; }
    Self& operator*=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self * v); return self; }
    Self& operator/=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self / v); return self; }
    Self& operator%=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self % v); return self; }
    Self& operator&=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self & v); return self; }
    Self& operator|=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self | v); return self; }
    Self& operator^=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self ^ v); return self; }
    Self& operator>>=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self >> v); return self; }
    Self& operator<<=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self << v); return self; }
};

template<class Self, class Other>
class IntegerBaseWith64: public ConstIntegerBaseWith64<Self, Other> {
private:
    inline static void assign(Self& self, uint32_t v) {
        IntegerBaseTrait<Self>::Dispatch::template assign<Self, uint32_t>(self, v);
    }

    inline static void assign(Self& self, uint64_t v) {
        IntegerBaseTrait<Self>::Dispatch::template assign<Self, uint64_t>(self, v);
    }
public:
    Self& operator+=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self + v); return self; }
    Self& operator-=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self - v); return self; }
    Self& operator*=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self * v); return self; }
    Self& operator/=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self / v); return self; }
    Self& operator%=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self % v); return self; }
    Self& operator&=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self & v); return self; }
    Self& operator|=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self | v); return self; }
    Self& operator^=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self ^ v); return self; }
    Self& operator>>=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self >> v); return self; }
    Self& operator<<=(const Other& v) { auto& self = static_cast<Self&>(*this); assign(self, self << v); return self; }
};

template<class Self>
class IntegerBaseWithSelf: public ConstIntegerBaseWithSelf<Self> {
private:
    typedef typename ConstIntegerBaseTrait<Self>::Dispatch::SelfMaxBit SelfMaxBit;
    inline static SelfMaxBit cast_for_self_op(const Self& self) {
        return ConstIntegerBaseTrait<Self>::Dispatch::template cast_for_op<Self, SelfMaxBit>(self);
    }

    inline static void assign(Self& self, uint32_t v) {
        IntegerBaseTrait<Self>::Dispatch::template assign<Self, uint32_t>(self, v);
    }

    inline static void assign(Self& self, uint64_t v) {
        IntegerBaseTrait<Self>::Dispatch::template assign<Self, uint64_t>(self, v);
    }
public:
    Self& operator+=(const Self& v)     { auto& self = static_cast<Self&>(*this); assign(self, self + v); return self;  }
    Self& operator-=(const Self& v)     { auto& self = static_cast<Self&>(*this); assign(self, self - v); return self;  }
    Self& operator*=(const Self& v)     { auto& self = static_cast<Self&>(*this); assign(self, self * v); return self;  }
    Self& operator/=(const Self& v)     { auto& self = static_cast<Self&>(*this); assign(self, self / v); return self;  }
    Self& operator%=(const Self& v)     { auto& self = static_cast<Self&>(*this); assign(self, self % v); return self;  }

    Self& operator++() { auto& self = static_cast<Self&>(*this); assign(self, self + 1u); return self; }
    Self& operator--() { auto& self = static_cast<Self&>(*this); assign(self, self - 1u); return self; }
    SelfMaxBit operator++(int) { auto& self = static_cast<Self&>(*this); auto tmp = cast_for_self_op(self); assign(self, self + 1u); return tmp; }
    SelfMaxBit operator--(int) { auto& self = static_cast<Self&>(*this); auto tmp = cast_for_self_op(self); assign(self, self - 1u); return tmp; }

    Self& operator&=(const Self& v)     { auto& self = static_cast<Self&>(*this); assign(self, self & v); return self;  }
    Self& operator|=(const Self& v)     { auto& self = static_cast<Self&>(*this); assign(self, self | v); return self;  }
    Self& operator^=(const Self& v)     { auto& self = static_cast<Self&>(*this); assign(self, self ^ v); return self;  }
    Self& operator>>=(const Self& v)    { auto& self = static_cast<Self&>(*this); assign(self, self >> v); return self; }
    Self& operator<<=(const Self& v)    { auto& self = static_cast<Self&>(*this); assign(self, self << v); return self; }
};


// NB: This is a hack needed to get the lookup rules for the operator memebers
// play nice. The basic issue is that different overloads in different base
// classes will be considered ambigious, so we need to import them all into
// the same base class.

template<typename...>
class IntegerBaseCombiner;

template<typename T, typename... Ts>
class IntegerBaseCombiner<T, Ts...>: public T, public IntegerBaseCombiner<Ts...> {
public:
    using T::operator +=;
    using T::operator -=;
    using T::operator *=;
    using T::operator /=;
    using T::operator %=;
    using T::operator &=;
    using T::operator |=;
    using T::operator ^=;
    using T::operator >>=;
    using T::operator <<=;

    using IntegerBaseCombiner<Ts...>::operator +=;
    using IntegerBaseCombiner<Ts...>::operator -=;
    using IntegerBaseCombiner<Ts...>::operator *=;
    using IntegerBaseCombiner<Ts...>::operator /=;
    using IntegerBaseCombiner<Ts...>::operator %=;
    using IntegerBaseCombiner<Ts...>::operator &=;
    using IntegerBaseCombiner<Ts...>::operator |=;
    using IntegerBaseCombiner<Ts...>::operator ^=;
    using IntegerBaseCombiner<Ts...>::operator >>=;
    using IntegerBaseCombiner<Ts...>::operator <<=;
};

template<typename T>
class IntegerBaseCombiner<T>: public T {
public:
    using T::operator +=;
    using T::operator -=;
    using T::operator *=;
    using T::operator /=;
    using T::operator %=;
    using T::operator &=;
    using T::operator |=;
    using T::operator ^=;
    using T::operator >>=;
    using T::operator <<=;
};

template<class Self>
using IntegerBase = IntegerBaseCombiner<
    IntegerBaseWithSelf<Self>,
    IntegerBaseWith32<Self, unsigned char>,
    IntegerBaseWith32<Self, char>,
    IntegerBaseWith32<Self, signed char>,
    IntegerBaseWith32<Self, unsigned short int>,
    IntegerBaseWith32<Self, signed short int>,
    IntegerBaseWith32<Self, unsigned int>,
    IntegerBaseWith32<Self, signed int>,
    IntegerBaseWith64<Self, unsigned long int>,
    IntegerBaseWith64<Self, signed long int>,
    IntegerBaseWith64<Self, unsigned long long int>,
    IntegerBaseWith64<Self, signed long long int>
>;

}

