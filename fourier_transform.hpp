#ifndef fourier_transform_h
#define fourier_transform_h
#include<vector>
#include<bit>
#include<algorithm>
#include<numbers>
#include<cmath>
#include<iterator>
#include<array>
#include<numeric>
#include<thread>
#include<future>
#include<type_traits>
namespace  yyk_fft {
constexpr unsigned int  yyk_N=10;
constexpr unsigned int yyk_thread_max=16;
template<typename C>
static std::array<std::array<C, (1<<yyk_N)>, yyk_thread_max> yyk_buffer {};
namespace  information {
template<typename U>
const static auto hc=static_cast<U>(std::countr_zero(std::bit_floor(std::min(std::max((unsigned int)2,std::thread::hardware_concurrency()),yyk_thread_max))));
static std::array<std::future<void>,yyk_thread_max> v{};
}
template<typename T,typename T1>
inline T byteswap(T x,T1 l) noexcept
{
    T y=0;
    for(T1 i=0;i<l;i++)
    {
        y=(y<<1)|(x&1);
        x>>=1;
    }
    return y;
}
template<typename InpuIt1,typename InputIt2,typename InputIt3,typename OutputIt,typename Tp>
inline OutputIt transform(InpuIt1 i1,InpuIt1 l1,InputIt2 i2,InputIt3 i3,OutputIt o,Tp tp) noexcept
{
    for(;i1!=l1;i1++)
        *o++=tp(*i1,*i2++,*i3++);
    return o;
}
template<bool parallel,bool Inverse,bool R,typename InputIt>
inline auto CooleyTukey(InputIt beg,InputIt end,const unsigned int N=0) noexcept
{
    using namespace information;
    using D=typename std::iterator_traits<InputIt>::difference_type;
    using C=typename std::iterator_traits<InputIt>::value_type;
    using U=std::conditional_t<sizeof(D)==8,unsigned  long long int ,unsigned  long int>;
    using T=typename C::value_type;
    constexpr T cs=Inverse?-2*std::numbers::pi_v<T>:2*std::numbers::pi_v<T>;
    auto& omiga=yyk_buffer<C>[N];
    auto u_len=static_cast<U>(end-beg);
    auto n=static_cast<U>(std::countr_zero(u_len));
    if constexpr(R){
        constexpr U L=1<<yyk_N;
        if(n<=yyk_N)
        {
            auto &rev=yyk_buffer<U>[N];
            for(U i=0;i<(1<<n);i++)
                rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (n - 1));
            for(U i=0;i<(1<<n);i++)
                if(i<rev[i]) std::iter_swap(beg+i, beg+rev[i]);
        }else{
            if constexpr(parallel)
            {
                auto &rev=yyk_buffer<U>;
                auto j=std::min(n-yyk_N,hc<U>);
                auto Rev=[&](auto i0,auto i1,auto d)
                {
                    for(auto i=i0;i<i1;i++){
                        for(U k=0;k<L;k++)
                            rev[d][k]=byteswap(L*i+k, n);
                        for(U k=0;k<L;k++) if((L*i+k)<rev[d][k])  std::iter_swap(beg+L*i+k, beg+rev[d][k]);
                    }
                };
                for(U i=0;i<(1<<j);i++)
                    v[i]=std::async(std::launch::async, Rev, i*(1<<(n-yyk_N-j)),(i+1)*(1<<(n-yyk_N-j)),i);
                std::for_each(v.begin(), v.begin()+(1<<j), [](auto &f){f.get();});
            }else{
                auto &rev=yyk_buffer<U>[N];
                auto j=n-yyk_N;
                for(U i=0;i<(1<<j);i++)
                {
                    for(U k=0;k<L;k++)
                        rev[k]=byteswap(i*L+k, n);
                    for(U k=0;k<L;k++)
                        if((L*i+k)<rev[k]) std::iter_swap(beg+L*i+k, beg+rev[k]);
                }
            }
        }
    }
    auto f=[](U i0,U i1,U t,U len,auto lbeg,auto rbeg,auto omiga_begin){
        for(U i=i0;i<i1;i++)
        {
             transform(lbeg, lbeg+len, rbeg, omiga_begin, rbeg, [](auto &l,auto &r,auto &m){
                auto _(l);
                auto __(m*r);
                l=_+__;
                return _-__;
            });
            rbeg=rbeg+(1<<(t+1));
            lbeg=lbeg+(1<<(t+1));
        }
    };
    auto g=[&](U t1){
        omiga[0]=C(1,0);
        for(U t=0;t<t1;t++)
        {
            C q(std::cos(cs/(1<<(t+1))),std::sin(cs/(1<<(t+1))));
            std::transform(omiga.begin(),std::prev(omiga.begin()+(1<<t)), std::next(omiga.begin()), [&q](auto &l){return q*l;});
            auto lbeg=beg;
            auto rbeg=lbeg+(1<<t);
            if constexpr(parallel)
            {
                auto l=std::min(n-t-1,hc<U>);
            for(U i=0;i<(1<<l);i++)
                v[i]=std::async(std::launch::async, f, i*(1<<(n-t-1-l)),(i+1)*(1<<(n-t-1-l)),t,1<<t,lbeg+i*(1<<(n-l)),rbeg+i*(1<<(n-l)),omiga.begin());
                std::for_each(v.begin(), v.begin()+(1<<l), [](auto &g){g.get();});
            }
            else {
                f(0,1<<(n-t-1),t,1<<t,lbeg,rbeg,omiga.begin());
            }
        }
    };
    if (n<=yyk_N+1)
    {
        g(n);
    }else{
        g(yyk_N+1);
        if constexpr(parallel)
        {
            auto &Omiga=yyk_buffer<C>;
            auto pf=[&](U j0,U j1,U i,U t,auto q)
            {
                for(U j=j0;j<j1;j++)
                {
                    auto index=j*(1<<yyk_N);
                    Omiga[i][0]=C(std::cos(cs*index/(1<<(t+1))),std::sin(cs*index/(1<<(t+1))));
                    std::transform(Omiga[i].begin(),std::prev(Omiga[i].end()), std::next(Omiga[i].begin()), [&q](auto &l){return q*l;});
                    auto lbeg=beg+index;
                    auto rbeg=lbeg+(1<<t);
                    f(0,1<<(n-t-1),t,1<<yyk_N,lbeg,rbeg,Omiga[i].begin());
                }
            };
            for(U t=yyk_N+1;t<n;t++)
            {
                C q(std::cos(cs/(1<<(t+1))),std::sin(cs/(1<<(t+1))));
                auto w=t-yyk_N;
                auto l=std::min(hc<U>,w);
                for(U i=0;i<(1<<l);i++)
                    v[i]=std::async(std::launch::async, pf, i*(1<<(w-l)),(i+1)*(1<<(w-l)),i,t,q);
                for(U i=0;i<(1<<l);i++) v[i].get();
            }
        }else{
            for(U t=yyk_N+1;t<n;t++)
            {
                C q(std::cos(cs/(1<<(t+1))),std::sin(cs/(1<<(t+1))));
                auto w=t-yyk_N;
                for(U j=0;j<(1<<w);j++)
                {
                    auto index=j*(1<<yyk_N);
                    omiga[0]=C(std::cos(cs*index/(1<<(t+1))),std::sin(cs*index/(1<<(t+1))));
                    std::transform(omiga.begin(),std::prev(omiga.end()), std::next(omiga.begin()), [&q](auto &l){return q*l;});
                    auto lbeg=beg+index;
                    auto rbeg=lbeg+(1<<t);
                    f(0,1<<(n-t-1),t,1<<yyk_N,lbeg,rbeg,omiga.begin());
                }
            }
        }
    }
    if constexpr(Inverse)
        std::for_each(beg, end, [&](auto &z){z/=u_len;});
    return beg;
}
template<typename T,typename F>
inline void naive_parallel_for(F f,T l,T r) noexcept
{
    auto len=r-l;
    auto th=information::hc<T>;
    auto d=len/th;
    for(T i=0;i<th-1;i++)
        information::v[i]=std::async(std::launch::async, f, l+d*i,l+d*(i+1));
    information::v[th-1]=std::async(std::launch::async, f, l+d*(th-1),r);
    for(T i=0;i<th;i++) information::v[i].get();
}
template<bool parallel,bool Inverse,typename InputIt>
inline InputIt Bluestein(InputIt beg,InputIt end) noexcept
{
    using C=typename std::iterator_traits<InputIt>::value_type;
    using T=typename C::value_type;
    constexpr T cs=Inverse?-std::numbers::pi_v<T>:std::numbers::pi_v<T>;
    using D=typename std::iterator_traits<InputIt>::difference_type;
    using U=std::conditional_t<sizeof(D)==8,unsigned  long long int ,unsigned  long int>;
    auto len=static_cast<U>(end-beg);
    T r=len&1?-1:1;
    auto u_len=std::bit_ceil<U>(2*len);
    auto n=static_cast<U>(std::countr_zero(u_len));
    std::vector<C> x(u_len);
    std::vector<C> y(u_len);
    std::vector<U> rev(len);
     for(U l=0;l<len;l++)
         rev[l]=(rev[l>>1]>>1)|((l&1)<<(n-1));
    auto calc_xy=[&](U l,U k){
    for(U i=l;i<k;i++)
    {
        auto alpha=i*i*cs/len;
        auto x0=std::cos(alpha);
        auto y0=std::sin(alpha);
        x[i]=C(x0,y0)*(*(beg+i));
        y[i]=C(r*x0,-r*y0);
    }
    };
    if constexpr(parallel)
        naive_parallel_for<U>(calc_xy, 0, len);
    else
        calc_xy(0,len);
    y[len]=C(1,0);
    std::reverse_copy(y.begin()+1, y.begin()+len, y.begin()+len+1);
    if constexpr(parallel)
    {
        auto f=[&](){
            for(U l=0;l<len;l++) if(l<rev[l]) {
                std::swap(x[l],x[rev[l]]);
                std::swap(y[l],y[rev[l]]);
            }
        };
        auto g=[&](){
            for(U l=len;l<2*len;l++) {
                auto b=(rev[l>>1]>>1)|((l&1)<<(n-1));
                if(l<b) std::swap(y[l], y[b]);
            }
        };
        auto t0=std::async(std::launch::async, f);
        auto t1=std::async(std::launch::async, g);
        t0.get();
        t1.get();
    }else{
        for(U l=0;l<len;l++) if(l<rev[l]) {
            std::swap(x[l],x[rev[l]]);
            std::swap(y[l],y[rev[l]]);
        }
        for(U l=len;l<2*len;l++) {
            auto b=(rev[l>>1]>>1)|((l&1)<<(n-1));
            if(l<b) std::swap(y[l], y[b]);
        }
    }
    CooleyTukey<parallel, false, false>(x.begin(), x.end(),0);
    CooleyTukey<parallel, false, false>(y.begin(), y.end(),1);
    std::transform(x.begin(), x.end(), y.begin(), x.begin(), [](auto &z1,auto &z2){return z1*z2;});
    CooleyTukey<parallel, true, true>(x.begin(), x.end());
    auto calc_f=[&](U l,U k){
    for(U i=l;i<k;i++)
    {
        auto alpha=i*i*cs/len;
        *(beg+i)=x[i+len]*C(std::cos(alpha),std::sin(alpha));
        if constexpr(Inverse)
            *(beg+i)/=len;
    }
    };
    if constexpr(parallel)
        naive_parallel_for<U>(calc_f, 0, len);
    else calc_f(0,len);
    return beg;
}
template<bool Inverse,typename InputIt>
inline InputIt dft(InputIt beg,InputIt end,const unsigned int n=0) noexcept
{
    using C=typename std::iterator_traits<InputIt>::value_type;
    using D=typename std::iterator_traits<InputIt>::difference_type;
    using T=typename C::value_type;
    constexpr T cs=Inverse?-2*std::numbers::pi_v<T>:2*std::numbers::pi_v<T>;
    auto len=end-beg;
    auto &m=yyk_buffer<C>[n];
    m[len]=std::accumulate(beg, end, C{});
    m[0]=C(1,0);
    for(D i=1;i<len;i++)
    {
        auto q(C(std::cos(cs*i/len),std::sin(cs*i/len)));
        std::transform(m.begin(), std::prev(m.begin()+len), std::next(m.begin()), [&](auto &l){return q*l;});
        m[len+i]=std::inner_product(beg, end, m.begin(), C{});
    }
    std::copy(m.begin()+len, m.begin()+2*len, beg);
    if constexpr(Inverse)
        std::for_each(beg, end, [&](auto &z){z/=len;});
    return beg;
}
template<int N,bool Inverse,typename InputIt>
inline constexpr InputIt special_dft(InputIt beg,const unsigned int n=0) noexcept
{
    static_assert(N==2||N==4||N==8, "N==2||N==4||N==8");
    using C=typename std::iterator_traits<InputIt>::value_type;
    using T=typename C::value_type;
    if constexpr(N==2)
    {
        auto x(*beg);
        auto y(*(beg+1));
        *beg=x+y;
        *(beg+1)=x-y;
    }
    else{
    if constexpr(N==4)
    {
        auto a(*beg);
        auto b(*(beg+1));
        auto c(*(beg+2));
        auto d(*(beg+3));
        *beg=a+b+c+d;
        *(beg+1)=a+C(0,Inverse?-1:1)*b+C(-1,0)*c+C(0,Inverse?1:-1)*d;
        *(beg+2)=a-b+c-d;
        *(beg+3)=a+C(0,Inverse?1:-1)*b+C(-1,0)*c+C(0,Inverse?-1:1)*d;
    }else{
      constexpr T  omiga[N][2] {{1., 0.}, {0.707107, 0.707107}, {0., 1.}, {-0.707107,
          0.707107}, {-1.,
          0.}, {-0.707107, -0.707107}, {0., -1.}, {0.707107, -0.707107}};
        std::array<C, 2*N> m;
        m[0]=C(1,0);
        m[N]=std::accumulate(beg, beg+N, C{});
        for(int i=1;i<N;i++)
        {
            C q(omiga[i][0],Inverse?-omiga[i][1]:omiga[i][1]);
            std::transform(m.begin(), std::prev(m.begin()+N), std::next(m.begin()), [&](auto &l){return q*l;});
            m[N+i]=std::inner_product(beg, beg+N, m.begin(), C{});
        }
        std::copy(m.begin()+N, m.begin()+2*N, beg);
    }
    }
    if constexpr(Inverse)
        std::for_each(beg,beg+N, [](auto &z){z/=N;});
    return beg;
}
template<bool parallel,bool Inverse,typename InputIt>
inline InputIt fourier_transform(InputIt beg,InputIt end) noexcept
{
    using T=typename std::iterator_traits<InputIt>::value_type;
    using D=typename std::iterator_traits<InputIt>::difference_type;
    using U=std::conditional_t<sizeof(D)==8,unsigned  long long int ,unsigned  long int>;
    auto len=end-beg;
    if(len<=1) return beg;
    static_assert(std::is_floating_point_v<typename T::value_type>, "Type Is Error");
    if(std::has_single_bit(static_cast<U>(len)))
    {
        switch (len) {
            case 2:
                return special_dft<2, Inverse>(beg);
            case 4:
                return special_dft<4, Inverse>(beg);
            case 8:
                return special_dft<8, Inverse>(beg);
            default:
                return CooleyTukey<parallel, Inverse, true>(beg, end);
        }
    }else{
        if(len<=40)
            return dft<Inverse>(beg, end);
        return Bluestein<parallel, Inverse>(beg, end);
    }
}
template<bool parallel=false,typename InputIt>
inline InputIt fast_fourier_transform(InputIt beg,InputIt end) noexcept
{
    return fourier_transform<parallel, false>(beg, end);
}
template<bool parallel=false,typename InputIt>
inline InputIt inverse_fast_fourier_transform(InputIt beg,InputIt end) noexcept
{
    return fourier_transform<parallel, true>(beg, end);
}
}
#endif /* fourier_transform_h */
