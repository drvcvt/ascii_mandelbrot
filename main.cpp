off")
        << "  aspect:" << (v.lock_aspect ? "lock" : "free")
        << "  frame:" << (int)ms << " ms"
        << "  prec(digits10):" << std::numeric_limits<Real>::digits10 << "\n"
        << "seq> " << std::flush;
}

// high-res PPM export
/**
 * Save a high-resolution PPM image of the current view.
 * Writes P6 binary format. Emits a message to std::cerr on I/O failure.
 */
static void save_ppm(const std::string &path, int W, int H, const View &v_in) {
    View v = v_in;
    v.ss = Cfg::SS_MIN;
    v.color = true;
    v.shade = true;

    std::ofstream f(path, std::ios::binary);
    if (!f) {
        std::cerr << "[save_ppm] failed to open '" << path << "' for writing\n";
        return;
    }

    f << "P6\n" << W << " " << H << "\n255\n";
    if (!f) {
        std::cerr << "[save_ppm] failed while writing header to '" << path << "'\n";
        return;
    }

    Real re_span = v.re_max - v.re_min;
    Real im_span = v.im_max - v.im_min;
    if (v.lock_aspect) {
        const Real px_aspect = Cfg::PX_ASPECT;
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
                Real boost = 1 - std::exp(-6.0L * std::max((Real)0.0, (Real)trap));
                t = std::clamp((t * 0.6L + boost * 0.5L), (Real)0, (Real)1);
            }
            Real br = v.shade ? shade_mu_gradient(re, im, v.max_it, dRe * Cfg::SHADE_DX_FACTOR, dIm * Cfg::SHADE_DX_FACTOR) : 1.0L;
            br = std::pow((double)br, Cfg::SHADE_GAMMA);
            Real h = std::fmod(t * Cfg::HUE_SCALE + Cfg::HUE_OFFSET, 1.0L);
            uint8_t R, G, B;
            hsv2rgb(h, Cfg::SATURATION, std::clamp(Cfg::VALUE_BASE + Cfg::VALUE_SCALE * br, 0.0L, 1.0L), R, G, B);
            f.put((char)R).put((char)G).put((char)B);
            if (!f) {
                std::cerr << "[save_ppm] write error while writing pixels to '" << path << "'\n";
                return;
            }
        }
    }
    f.close();
    if (!f) {
        std::cerr << "[save_ppm] close/error flushing '" << path << "'\n";
    }
}

// ---- main + event queue/step-anim ----
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    const int W = Cfg::TERM_W, H = Cfg::TERM_H;
    View view;
    bool step_anim = false; // y toggelt
    int step_delay_ms = Cfg::STEP_DELAY_MS_DEFAULT; // , / . ändern das

    auto apply_key = [&](char c) -> bool { // return true wenn state geändert/handled
        Real re_span = view.re_max - view.re_min;
        Real im_span = view.im_max - view.im_min;
        Real re_c = (view.re_min + view.re_max) * 0.5L;
        Real im_c = (view.im_min + view.im_max) * 0.5L;
        switch (c) {
        case '+': {
            re_span *= Cfg::ZOOM_FACTOR_IN;
            im_span *= Cfg::ZOOM_FACTOR_IN;
            view.re_min = re_c - re_span * 0.5L;
            view.re_max = re_c + re_span * 0.5L;
            view.im_min = im_c - im_span * 0.5L;
            view.im_max = im_c + im_span * 0.5L;
            return true; }
        case '-': {
            re_span *= Cfg::ZOOM_FACTOR_OUT;
            im_span *= Cfg::ZOOM_FACTOR_OUT;
            view.re_min = re_c - re_span * 0.5L;
            view.re_max = re_c + re_span * 0.5L;
            view.im_min = im_c - im_span * 0.5L;
            view.im_max = im_c + im_span * 0.5L;
            return true; }
        case 'a': {
            Real dx = -Cfg::PAN_STEP * re_span;
            view.re_min += dx; view.re_max += dx; return true; }
        case 'd': {
            Real dx =  Cfg::PAN_STEP * re_span;
            view.re_min += dx; view.re_max += dx; return true; }
        case 'w': {
            Real dy =  Cfg::PAN_STEP * im_span;
            view.im_min += dy; view.im_max += dy; return true; }
        case 's': {
            Real dy = -Cfg::PAN_STEP * im_span;
            view.im_min += dy; view.im_max += dy; return true; }
        case '[':
            view.max_it = std::max(Cfg::ITER_MIN, view.max_it - Cfg::ITER_STEP);
            view.auto_iters = false; return true;
        case ']':
            view.max_it = std::min(Cfg::ITER_MAX, view.max_it + Cfg::ITER_STEP);
            view.auto_iters = false; return true;
        case 't':
            view.ss = (view.ss % Cfg::SS_MAX) + 1; return true;
        case 'i':
            view.auto_iters = !view.auto_iters; return true;
        case 'c':
            view.color = !view.color; return true;
        case 'g':
            view.shade = !view.shade; return true;
        case 'o':
            view.orbit = !view.orbit; return true;
        case 'z':
            view.lock_aspect = !view.lock_aspect; return true;
        case 'p':
            save_ppm(Cfg::PPM_PATH, Cfg::PPM_W, Cfg::PPM_H, view); return false;
        case 'r':
            view = View{}; return true;
        case ',':
            step_delay_ms = std::min(Cfg::STEP_DELAY_MS_MAX, step_delay_ms + Cfg::STEP_DELAY_MS_DELTA); return false;
        case '.':
            step_delay_ms = std::max(Cfg::STEP_DELAY_MS_MIN, step_delay_ms - Cfg::STEP_DELAY_MS_DELTA); return false;
        case 'y':
            step_anim = !step_anim; return false;
        case ' ': case '\t': case '\r': case '\n':
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
        if (seq.find('q') != std::string::npos || seq.find('Q') != std::string::npos) {
            clear_screen();
            reset_color();
            return 0;
        }
        // queue & process: char für char
        std::deque<char> q;
        for (char c : seq) q.push_back(c);
        if (step_anim) {
            while (!q.empty()) {
                char c = q.front(); q.pop_front();
                bool changed = apply_key(c);
                if (changed) {
                    render(W, H, view);
                    if (step_delay_ms > 0)
                        std::this_thread::sleep_for(std::chrono::milliseconds(step_delay_ms));
                }
            }
        } else {
            bool any = false;
            while (!q.empty()) { any |= apply_key(q.front()); q.pop_front(); }
            // next loop renders end-state
        }
    }
}
