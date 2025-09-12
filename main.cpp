// mandel_deep_queue.cpp
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <thread>

using Real = long double;

static const std::string PALETTE =
    " .'`^\",:;Il!i~+_-?][}{1)(|\\/*tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";

// --- input: ganze zeile (ENTER nötig), wir verarbeiten jeden char ---
static std::string read_line() {
    std::string s;
    if (!std::getline(std::cin, s))
        return "q"; // EOF -> quit
    return s;
}

struct View {
    Real re_min = -2.5, re_max = 1.0;
    Real im_min = -1.2, im_max = 1.2;
    int max_it = 240;
    int ss = 1;              // supersampling: 1/2/3
    bool color = true;       // 24-bit ANSI
    bool shade = true;       // hill-shading
    bool orbit = false;      // orbit-trap mix
    bool lock_aspect = true; // fix terminal aspect
    bool auto_iters = true;  // iters nach zoom
};

static void clear_screen() { std::cout << "\x1b[2J\x1b[H"; }
static void reset_color() { std::cout << "\x1b[0m"; }

// HSV -> RGB (0..1)
static void hsv2rgb(Real h, Real s, Real v, uint8_t &R, uint8_t &G,
                    uint8_t &B) {
    h = std::fmod(std::max<Real>(0, h), (Real)1.0) * 6.0L;
    Real c = v * s;
    Real x = c * (1 - fabsl(fmodl(h, 2.0L) - 1));
    Real m = v - c;
    Real r = 0, g = 0, b = 0;
    int i = (int)h;
    switch (i) {
    case 0:
        r = c;
        g = x;
        b = 0;
        break;
    case 1:
        r = x;
        g = c;
        b = 0;
        break;
    case 2:
        r = 0;
        g = c;
        b = x;
        break;
    case 3:
        r = 0;
        g = x;
        b = c;
        break;
    case 4:
        r = x;
        g = 0;
        b = c;
        break;
    default:
        r = c;
        g = 0;
        b = x;
        break;
    }
    R = (uint8_t)std::clamp((int)((r + m) * 255 + 0.5), 0, 255);
    G = (uint8_t)std::clamp((int)((g + m) * 255 + 0.5), 0, 255);
    B = (uint8_t)std::clamp((int)((b + m) * 255 + 0.5), 0, 255);
}
static void set_rgb(uint8_t R, uint8_t G, uint8_t B) {
    std::cout << "\x1b[38;2;" << (int)R << ';' << (int)G << ';' << (int)B
              << 'm';
}

static int auto_iters(Real re_span) {
    Real s = std::log10((double)re_span);
    int base = 240, per_dec = 150;
    int it = base + (int)(-s * per_dec);
    return std::clamp(it, 120, 20000);
}

static inline bool inside_fast(Real x, Real y) {
    Real yy = y * y;
    Real xq = (x - 0.25L);
    Real q = xq * xq + yy;
    if (q * (q + xq) <= 0.25L * yy)
        return true; // cardioid
    if ((x + 1.0L) * (x + 1.0L) + yy <= 0.0625L)
        return true; // period-2 bulb
    return false;
}

struct Stats {
    int it = 0;
    Real zr = 0, zi = 0;
    bool inside = false;
    Real min_r = 1e9L;
    Real min_axis = 1e9L;
};

static inline Stats iterate_point(Real re, Real im, int max_it) {
    Stats s{};
    if (inside_fast(re, im)) {
        s.it = max_it;
        s.inside = true;
        return s;
    }
    Real zr = 0, zi = 0;
    int it = 0;
    for (; it < max_it; ++it) {
        Real r = std::hypotl(zr, zi);
        if (r < s.min_r)
            s.min_r = r;
        Real a = fabsl(zi);
        if (a < s.min_axis)
            s.min_axis = a;

        Real zr2 = zr * zr - zi * zi + re;
        Real zi2 = 2 * zr * zi + im;
        zr = zr2;
        zi = zi2;
        if (zr * zr + zi * zi > 4) {
            ++it;
            break;
        }
    }
    s.it = it;
    s.zr = zr;
    s.zi = zi;
    s.inside = (it >= max_it);
    return s;
}

static inline Real smooth_mu(const Stats &s) {
    if (s.inside)
        return 0;
    Real r2 = s.zr * s.zr + s.zi * s.zi;
    r2 = std::max<Real>(r2, 1.0L + 1e-12L);
    Real mu =
        s.it - (std::log(std::log(std::sqrt((double)r2))) / std::log(2.0L));
    return std::isfinite((double)mu) ? mu : (Real)s.it;
}

static inline Real shade_mu_gradient(Real re, Real im, int max_it, Real dx,
                                     Real dy) {
    Stats sx1 = iterate_point(re + dx, im, max_it);
    Stats sx2 = iterate_point(re - dx, im, max_it);
    Stats sy1 = iterate_point(re, im + dy, max_it);
    Stats sy2 = iterate_point(re, im - dy, max_it);
    Real mux = smooth_mu(sx1) - smooth_mu(sx2);
    Real muy = smooth_mu(sy1) - smooth_mu(sy2);
    Real nx = -mux, ny = -muy, nz = 1.0L;
    Real inv =
        1.0L / std::max<Real>(1e-12L, std::sqrt(nx * nx + ny * ny + nz * nz));
    nx *= inv;
    ny *= inv;
    nz *= inv;
    Real lx = -0.45L, ly = 0.55L, lz = 0.65L;
    Real linv = 1.0L / std::sqrt(lx * lx + ly * ly + lz * lz);
    lx *= linv;
    ly *= linv;
    lz *= linv;
    Real dot = nx * lx + ny * ly + nz * lz;
    return std::clamp(dot, 0.0L, 1.0L);
}

static void render(int W, int H, View &v) {
    if (v.auto_iters)
        v.max_it = auto_iters(v.re_max - v.re_min);

    clear_screen();
    auto t0 = std::chrono::high_resolution_clock::now();

    Real re_span = v.re_max - v.re_min;
    Real im_span = v.im_max - v.im_min;
    Real im_min_r = v.im_min, im_max_r = v.im_max;
    if (v.lock_aspect) {
        const Real char_aspect = 2.0L;
        Real want = re_span * ((Real)H / (Real)W) * char_aspect;
        Real ic = (v.im_min + v.im_max) * 0.5L;
        im_min_r = ic - want * 0.5L;
        im_max_r = ic + want * 0.5L;
        im_span = want;
    }

    Real dRe = re_span / (W - 1);
    Real dIm = im_span / (H - 1);

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            int samples = v.ss * v.ss;
            Real acc_t = 0, acc_b = 0;
            bool all_inside = true;
            Real re0 = v.re_min + x * dRe;
            Real im0 = im_max_r - y * dIm;

            for (int sy = 0; sy < v.ss; ++sy) {
                for (int sx = 0; sx < v.ss; ++sx) {
                    Real re = re0 + ((sx + 0.5L) / v.ss - 0.5L) * dRe;
                    Real im = im0 - ((sy + 0.5L) / v.ss - 0.5L) * dIm;

                    Stats s = iterate_point(re, im, v.max_it);
                    if (!s.inside)
                        all_inside = false;

                    Real mu = smooth_mu(s);
                    Real t = std::clamp(mu / v.max_it, (Real)0, (Real)1);

                    if (v.orbit) {
                        Real trap = std::min(s.min_r, s.min_axis);
                        Real boost = 1 - std::exp(-6.0L * std::max((Real)0.0,
                                                                   (Real)trap));
                        t = std::clamp((t * 0.6L + boost * 0.5L), (Real)0,
                                       (Real)1);
                    }

                    Real b = 1.0L;
                    if (v.shade) {
                        Real br = shade_mu_gradient(re, im, v.max_it,
                                                    dRe * 0.7L, dIm * 0.7L);
                        b = std::pow((double)br, 0.8);
                    }

                    acc_t += t;
                    acc_b += b;
                }
            }

            Real t_avg = acc_t / samples;
            Real b_avg = std::clamp(acc_b / samples, (Real)0, (Real)1);

            if (v.color) {
                Real h = std::fmod(t_avg * 0.85L + 0.02L, 1.0L);
                Real sS = 0.9L;
                Real vV = std::clamp(0.15L + 0.95L * b_avg, 0.0L, 1.0L);
                uint8_t R, G, B;
                hsv2rgb(h, sS, vV, R, G, B);
                set_rgb(R, G, B);
                std::cout << "█";
            } else {
                int idx = (int)std::clamp(
                    (int)(1 + (PALETTE.size() - 2) * (t_avg * b_avg)), 1,
                    (int)PALETTE.size() - 1);
                if (all_inside)
                    idx = 0;
                std::cout << PALETTE[idx];
            }
        }
        if (v.color)
            reset_color();
        std::cout << '\n';
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout
        << "\n[controls]  type sequences, ENTER (e.g. ddd++a-)\n"
        << "  +/- zoom   wasd pan   [ ] iters   t supersample   i auto-iters   "
           "c color   g shade   o orbit   z aspect   p save PPM   y step-anim  "
           " ,/. anim speed   r reset   q quit\n"
        << "re:[" << (double)v.re_min << ", " << (double)v.re_max << "]  "
        << "im:[" << (double)v.im_min << ", " << (double)v.im_max << "]  "
        << "iters:" << v.max_it << "  auto:" << (v.auto_iters ? "on" : "off")
        << "  ss:" << v.ss << "  color:" << (v.color ? "on" : "off")
        << "  shade:" << (v.shade ? "on" : "off")
        << "  orbit:" << (v.orbit ? "on" : "off")
        << "  aspect:" << (v.lock_aspect ? "lock" : "free")
        << "  frame:" << (int)ms << " ms"
        << "  prec(digits10):" << std::numeric_limits<Real>::digits10 << "\n"
        << "seq> " << std::flush;
}

// high-res PPM export
static void save_ppm(const std::string &path, int W, int H, const View &v_in) {
    View v = v_in;
    v.ss = 1;
    v.color = true;
    v.shade = true;
    std::ofstream f(path, std::ios::binary);
    if (!f)
        return;
    f << "P6\n" << W << " " << H << "\n255\n";

    Real re_span = v.re_max - v.re_min;
    Real im_span = v.im_max - v.im_min;
    if (v.lock_aspect) {
        const Real px_aspect = 1.0L;
        Real want = re_span * ((Real)H / (Real)W) * px_aspect;
        Real ic = (v.im_min + v.im_max) * 0.5L;
        v.im_min = ic - want * 0.5L;
        v.im_max = ic + want * 0.5L;
        im_span = want;
    }
    Real dRe = re_span / (W - 1);
    Real dIm = im_span / (H - 1);

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Real re = v.re_min + x * dRe;
            Real im = v.im_max - y * dIm;
            Stats s = iterate_point(re, im, v.max_it);
            Real t = std::clamp(smooth_mu(s) / v.max_it, (Real)0, (Real)1);
            if (v.orbit) {
                Real trap = std::min(s.min_r, s.min_axis);
                Real boost =
                    1 - std::exp(-6.0L * std::max((Real)0.0, (Real)trap));
                t = std::clamp((t * 0.6L + boost * 0.5L), (Real)0, (Real)1);
            }
            Real br = v.shade ? shade_mu_gradient(re, im, v.max_it, dRe * 0.7L,
                                                  dIm * 0.7L)
                              : 1.0L;
            br = std::pow((double)br, 0.8);
            Real h = std::fmod(t * 0.85L + 0.02L, 1.0L);
            uint8_t R, G, B;
            hsv2rgb(h, 0.9L, std::clamp(0.15L + 0.95L * br, 0.0L, 1.0L), R, G,
                    B);
            f.put((char)R).put((char)G).put((char)B);
        }
    }
    f.close();
}

// ---- main + event queue/step-anim ----
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    const int W = 120, H = 40;
    View view;

    bool step_anim = false; // y toggelt
    int step_delay_ms = 40; // , / . ändern das

    auto apply_key =
        [&](char c) -> bool { // return true wenn state geändert/handled
        Real re_span = view.re_max - view.re_min;
        Real im_span = view.im_max - view.im_min;
        Real re_c = (view.re_min + view.re_max) * 0.5L;
        Real im_c = (view.im_min + view.im_max) * 0.5L;

        switch (c) {
        case '+':
            re_span *= 0.5L;
            im_span *= 0.5L;
            view.re_min = re_c - re_span * 0.5L;
            view.re_max = re_c + re_span * 0.5L;
            view.im_min = im_c - im_span * 0.5L;
            view.im_max = im_c + im_span * 0.5L;
            return true;
        case '-':
            re_span *= 2.0L;
            im_span *= 2.0L;
            view.re_min = re_c - re_span * 0.5L;
            view.re_max = re_c + re_span * 0.5L;
            view.im_min = im_c - im_span * 0.5L;
            view.im_max = im_c + im_span * 0.5L;
            return true;
        case 'a': {
            Real dx = -0.15L * re_span;
            view.re_min += dx;
            view.re_max += dx;
            return true;
        }
        case 'd': {
            Real dx = 0.15L * re_span;
            view.re_min += dx;
            view.re_max += dx;
            return true;
        }
        case 'w': {
            Real dy = 0.15L * im_span;
            view.im_min += dy;
            view.im_max += dy;
            return true;
        }
        case 's': {
            Real dy = -0.15L * im_span;
            view.im_min += dy;
            view.im_max += dy;
            return true;
        }
        case '[':
            view.max_it = std::max(10, view.max_it - 100);
            view.auto_iters = false;
            return true;
        case ']':
            view.max_it = std::min(20000, view.max_it + 100);
            view.auto_iters = false;
            return true;
        case 't':
            view.ss = (view.ss % 3) + 1;
            return true;
        case 'i':
            view.auto_iters = !view.auto_iters;
            return true;
        case 'c':
            view.color = !view.color;
            return true;
        case 'g':
            view.shade = !view.shade;
            return true;
        case 'o':
            view.orbit = !view.orbit;
            return true;
        case 'z':
            view.lock_aspect = !view.lock_aspect;
            return true;
        case 'p':
            save_ppm("mandel.ppm", 1024, 768, view);
            return false;
        case 'r':
            view = View{};
            return true;
        case ',':
            step_delay_ms = std::min(400, step_delay_ms + 20);
            return false;
        case '.':
            step_delay_ms = std::max(0, step_delay_ms - 20);
            return false;
        case 'y':
            step_anim = !step_anim;
            return false;
        case ' ':
        case '\t':
        case '\r':
        case '\n':
            return false;
        default:
            return false;
        }
    };

    while (true) {
        render(W, H, view);
        std::string seq = read_line();
        if (seq.empty())
            continue;

        // simple quit, anywhere in sequence
        if (seq.find('q') != std::string::npos ||
            seq.find('Q') != std::string::npos) {
            clear_screen();
            reset_color();
            return 0;
        }

        // queue & process: char für char
        std::deque<char> q;
        for (char c : seq)
            q.push_back(c);

        if (step_anim) {
            while (!q.empty()) {
                char c = q.front();
                q.pop_front();
                bool changed = apply_key(c);
                if (changed) {
                    render(W, H, view);
                    if (step_delay_ms > 0)
                        std::this_thread::sleep_for(
                            std::chrono::milliseconds(step_delay_ms));
                }
            }
        } else {
            bool any = false;
            while (!q.empty()) {
                any |= apply_key(q.front());
                q.pop_front();
            }
            // nächster loop rendert dann den Endzustand
        }
    }
}
