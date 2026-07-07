// VortexStreet — 2D lattice Boltzmann (D2Q9, BGK) channel flow
// Kármán vortex street behind a cylinder, passive dye advection,
// mouse-drawn obstacles, dye / speed / vorticity / density views.
// C++17 / Qt6 Widgets, no other dependencies.

#include <QApplication>
#include <QWidget>
#include <QPainter>
#include <QTimer>
#include <QElapsedTimer>
#include <QMouseEvent>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QComboBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <QLabel>
#include <QGroupBox>
#include <QStyleFactory>

#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <thread>
#include <vector>

static inline float clampf(float v, float lo, float hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

// D2Q9 lattice
static constexpr int Q = 9;
static constexpr int EX[Q] = { 0, 1, 0, -1,  0, 1, -1, -1,  1 };
static constexpr int EY[Q] = { 0, 0, 1,  0, -1, 1,  1, -1, -1 };
static constexpr int OPP[Q] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
static constexpr float W[Q] = {
    4.f / 9.f,
    1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f,
    1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f
};

enum class View { Dye, Speed, Vorticity, Density };
enum class Collision { BGK, MRT };

struct SimParams {
    Collision collision = Collision::MRT;
    float omega        = 1.92f;  // relaxation, viscosity nu = (1/omega - 0.5)/3
    float u0           = 0.09f;  // inflow speed, lattice units
    int   stepsPerFrame = 12;
    float cylinderR    = 0.f;    // set on reset from grid size
    View  view         = View::Dye;
    bool  injectDye    = true;
    bool  paused       = false;
};

class LbmCanvas : public QWidget
{
    Q_OBJECT
public:
    explicit LbmCanvas(QWidget *parent = nullptr)
    {
        setMinimumSize(640, 340);
        setGridPreset(1);
        connect(&m_timer, &QTimer::timeout, this, &LbmCanvas::step);
        m_timer.start(16);
    }

    SimParams params;

    void setGridPreset(int idx)
    {
        static const int dims[3][2] = { {400, 200}, {600, 300}, {900, 450} };
        m_nx = dims[idx][0];
        m_ny = dims[idx][1];
        reset(true);
    }

    void reset(bool rebuildSolids)
    {
        const size_t total = size_t(m_nx) * m_ny;
        for (int q = 0; q < Q; ++q) {
            m_f[q].assign(total, 0.f);
            m_fp[q].assign(total, 0.f);
        }
        m_u.assign(total, 0.f);
        m_v.assign(total, 0.f);
        m_rho.assign(total, 1.f);
        m_dye.assign(total, 0.f);
        m_dyeTmp.assign(total, 0.f);
        m_image = QImage(m_nx, m_ny, QImage::Format_RGB32);

        if (rebuildSolids) {
            m_solid.assign(total, 0);
            // channel walls
            for (int x = 0; x < m_nx; ++x) {
                m_solid[idx(x, 0)] = 1;
                m_solid[idx(x, m_ny - 1)] = 1;
            }
            // cylinder, slightly off-center to trigger shedding
            params.cylinderR = m_ny * 0.08f;
            addDisk(m_nx / 4, m_ny / 2 + 2, params.cylinderR, true);
        }

        // impulsive start: equilibrium at (u0, 0)
        const float u0 = params.u0;
        for (int y = 0; y < m_ny; ++y)
            for (int x = 0; x < m_nx; ++x) {
                size_t i = idx(x, y);
                for (int q = 0; q < Q; ++q)
                    m_f[q][i] = feq(q, 1.f, m_solid[i] ? 0.f : u0, 0.f);
            }
    }

    void clearObstacles()
    {
        const size_t total = size_t(m_nx) * m_ny;
        for (size_t i = 0; i < total; ++i) m_solid[i] = 0;
        for (int x = 0; x < m_nx; ++x) {
            m_solid[idx(x, 0)] = 1;
            m_solid[idx(x, m_ny - 1)] = 1;
        }
    }

    float msPerStep() const { return m_msPerStep; }
    float reynolds() const
    {
        float nu = (1.f / params.omega - 0.5f) / 3.f;
        return params.u0 * (2.f * params.cylinderR) / std::max(nu, 1e-6f);
    }

signals:
    void statsChanged();

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter p(this);
        p.fillRect(rect(), QColor(6, 6, 10));
        const float ar = float(m_nx) / m_ny;
        int w = width(), h = int(w / ar);
        if (h > height()) { h = height(); w = int(h * ar); }
        m_view = QRect((width() - w) / 2, (height() - h) / 2, w, h);
        p.setRenderHint(QPainter::SmoothPixmapTransform);
        p.drawImage(m_view, m_image);
    }

    void mousePressEvent(QMouseEvent *e) override { drawAt(e->pos(), e->buttons()); }
    void mouseMoveEvent(QMouseEvent *e) override  { drawAt(e->pos(), e->buttons()); }

private slots:
    void step()
    {
        if (!params.paused) {
            QElapsedTimer t; t.start();
            for (int s = 0; s < params.stepsPerFrame; ++s)
                lbmStep();
            m_msPerStep = float(t.nsecsElapsed()) * 1e-6f / params.stepsPerFrame;
            advectDye(float(params.stepsPerFrame));
        }
        composeImage();
        update();
        if (const char *df = getenv("DUMP_FRAMES")) {
            static int left = atoi(df);
            if (--left <= 0) { m_image.save("dump.png"); QApplication::quit(); }
        }
        if (!m_statClock.isValid() || m_statClock.hasExpired(400)) {
            emit statsChanged();
            m_statClock.restart();
        }
    }

private:
    size_t idx(int x, int y) const { return size_t(y) * m_nx + x; }

    static int threadCount()
    {
        return std::max(2u, std::thread::hardware_concurrency());
    }

    template <class F> void parallelRows(int y0, int y1, F &&fn) const
    {
        const int T = threadCount();
        std::vector<std::thread> th;
        for (int t = 0; t < T; ++t)
            th.emplace_back([&, t] { for (int y = y0 + t; y < y1; y += T) fn(y); });
        for (auto &x : th) x.join();
    }

    static inline float feq(int q, float rho, float ux, float uy)
    {
        float eu = EX[q] * ux + EY[q] * uy;
        float u2 = ux * ux + uy * uy;
        return W[q] * rho * (1.f + 3.f * eu + 4.5f * eu * eu - 1.5f * u2);
    }

    void addDisk(int cx, int cy, float r, bool solid)
    {
        const int ri = int(r) + 1;
        for (int y = std::max(1, cy - ri); y <= std::min(m_ny - 2, cy + ri); ++y)
            for (int x = std::max(0, cx - ri); x <= std::min(m_nx - 1, cx + ri); ++x) {
                float dx = x - cx, dy = y - cy;
                if (dx * dx + dy * dy <= r * r) {
                    size_t i = idx(x, y);
                    bool was = m_solid[i];
                    m_solid[i] = solid ? 1 : 0;
                    if (was && !solid) {
                        // re-fluidized cell: fill with local equilibrium
                        for (int q = 0; q < Q; ++q)
                            m_f[q][i] = feq(q, 1.f, 0.f, 0.f);
                        m_dye[i] = 0.f;
                    }
                }
            }
    }

    void drawAt(QPoint pos, Qt::MouseButtons b)
    {
        if (!m_view.contains(pos) || !(b & (Qt::LeftButton | Qt::RightButton)))
            return;
        int cx = int((pos.x() - m_view.x()) / float(m_view.width()) * m_nx);
        int cy = int((pos.y() - m_view.y()) / float(m_view.height()) * m_ny);
        if (cx < 6) return; // keep the inlet clean
        addDisk(cx, cy, std::max(3.f, m_ny * 0.02f), b & Qt::LeftButton);
    }

    void lbmStep()
    {
        const float omega = params.omega;
        const float u0 = params.u0;
        const int nx = m_nx, ny = m_ny;

        // 1) collide (into m_fp), store macroscopic fields
        parallelRows(0, ny, [&](int y) {
            for (int x = 0; x < nx; ++x) {
                const size_t i = idx(x, y);
                if (m_solid[i]) {
                    for (int q = 0; q < Q; ++q) m_fp[q][i] = m_f[q][i];
                    m_u[i] = m_v[i] = 0.f;
                    continue;
                }
                float rho = 0.f, ux = 0.f, uy = 0.f;
                float fq[Q];
                for (int q = 0; q < Q; ++q) {
                    fq[q] = m_f[q][i];
                    rho += fq[q];
                    ux += fq[q] * EX[q];
                    uy += fq[q] * EY[q];
                }
                float inv = 1.f / rho;
                ux *= inv; uy *= inv;
                // clamp for stability at high Re
                float sp = std::sqrt(ux * ux + uy * uy);
                if (sp > 0.35f) { ux *= 0.35f / sp; uy *= 0.35f / sp; }
                m_rho[i] = rho; m_u[i] = ux; m_v[i] = uy;
                if (params.collision == Collision::BGK) {
                    for (int q = 0; q < Q; ++q)
                        m_fp[q][i] = fq[q] + omega * (feq(q, rho, ux, uy) - fq[q]);
                } else {
                    // MRT (Lallemand & Luo). M rows are orthogonal; norms 9,36,36,6,12,6,12,4,4.
                    // Moment order: rho, e, eps, jx, qx, jy, qy, pxx, pxy.
                    static constexpr float M[9][Q] = {
                        { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                        {-4,-1,-1,-1,-1, 2, 2, 2, 2 },
                        { 4,-2,-2,-2,-2, 1, 1, 1, 1 },
                        { 0, 1, 0,-1, 0, 1,-1,-1, 1 },
                        { 0,-2, 0, 2, 0, 1,-1,-1, 1 },
                        { 0, 0, 1, 0,-1, 1, 1,-1,-1 },
                        { 0, 0,-2, 0, 2, 1, 1,-1,-1 },
                        { 0, 1,-1, 1,-1, 0, 0, 0, 0 },
                        { 0, 0, 0, 0, 0, 1,-1, 1,-1 }
                    };
                    static constexpr float invNorm[9] = {
                        1.f/9, 1.f/36, 1.f/36, 1.f/6, 1.f/12, 1.f/6, 1.f/12, 1.f/4, 1.f/4
                    };
                    float m[9];
                    for (int k = 0; k < 9; ++k) {
                        float acc = 0.f;
                        for (int q = 0; q < Q; ++q) acc += M[k][q] * fq[q];
                        m[k] = acc;
                    }
                    const float jx = rho * ux, jy = rho * uy;
                    const float j2r = 3.f * (jx * jx + jy * jy) / rho;
                    const float meq[9] = {
                        rho, -2.f * rho + j2r, rho - j2r,
                        jx, -jx, jy, -jy,
                        (jx * jx - jy * jy) / rho, jx * jy / rho
                    };
                    // relaxation rates: conserved moments 0; s_nu = omega sets viscosity
                    const float S[9] = { 0.f, 1.64f, 1.54f, 0.f, 1.9f, 0.f, 1.9f, omega, omega };
                    for (int k = 1; k < 9; ++k)
                        m[k] -= S[k] * (m[k] - meq[k]);
                    for (int q = 0; q < Q; ++q) {
                        float acc = 0.f;
                        for (int k = 0; k < 9; ++k)
                            acc += M[k][q] * invNorm[k] * m[k];
                        m_fp[q][i] = acc;
                    }
                }
            }
        });

        // 2) stream (pull) with half-way bounce-back on solids
        parallelRows(0, ny, [&](int y) {
            for (int x = 0; x < nx; ++x) {
                const size_t i = idx(x, y);
                if (m_solid[i]) continue;
                for (int q = 0; q < Q; ++q) {
                    int sx = x - EX[q], sy = y - EY[q];
                    if (sx < 0 || sx >= nx) { // handled by inlet/outlet below
                        m_f[q][i] = m_fp[q][i];
                        continue;
                    }
                    size_t s = idx(sx, sy);
                    m_f[q][i] = m_solid[s] ? m_fp[OPP[q]][i] : m_fp[q][s];
                }
            }
        });

        // 3) inlet: fixed equilibrium; outlet: zero-gradient copy
        parallelRows(1, ny - 1, [&](int y) {
            size_t iIn = idx(0, y);
            if (!m_solid[iIn])
                for (int q = 0; q < Q; ++q)
                    m_f[q][iIn] = feq(q, 1.f, u0, 0.f);
            size_t iOut = idx(nx - 1, y), iPrev = idx(nx - 2, y);
            if (!m_solid[iOut])
                for (int q = 0; q < Q; ++q)
                    m_f[q][iOut] = m_f[q][iPrev];
        });
    }

    void advectDye(float dt)
    {
        const int nx = m_nx, ny = m_ny;
        // inject stripes at the inlet
        if (params.injectDye) {
            const int stripe = std::max(6, ny / 24);
            for (int y = 1; y < ny - 1; ++y)
                if ((y / stripe) % 2 == 0)
                    for (int x = 1; x <= 3; ++x)
                        m_dye[idx(x, y)] = 1.f;
        }
        // semi-Lagrangian backtrace
        parallelRows(1, ny - 1, [&](int y) {
            for (int x = 1; x < nx - 1; ++x) {
                size_t i = idx(x, y);
                if (m_solid[i]) { m_dyeTmp[i] = 0.f; continue; }
                float bx = clampf(x - m_u[i] * dt, 0.5f, nx - 1.5f);
                float by = clampf(y - m_v[i] * dt, 0.5f, ny - 1.5f);
                int x0 = int(bx), y0 = int(by);
                float fx = bx - x0, fy = by - y0;
                float d00 = m_dye[idx(x0, y0)],     d10 = m_dye[idx(x0 + 1, y0)];
                float d01 = m_dye[idx(x0, y0 + 1)], d11 = m_dye[idx(x0 + 1, y0 + 1)];
                float d = (d00 * (1 - fx) + d10 * fx) * (1 - fy)
                        + (d01 * (1 - fx) + d11 * fx) * fy;
                m_dyeTmp[i] = d * 0.9995f;
            }
        });
        std::swap(m_dye, m_dyeTmp);
    }

    void composeImage()
    {
        const int nx = m_nx, ny = m_ny;
        const View view = params.view;
        const float invU = 1.f / std::max(params.u0 * 1.8f, 1e-4f);
        parallelRows(0, ny, [&](int y) {
            QRgb *out = reinterpret_cast<QRgb *>(m_image.scanLine(y));
            for (int x = 0; x < nx; ++x) {
                size_t i = idx(x, y);
                if (m_solid[i]) { out[x] = qRgb(120, 122, 128); continue; }
                float r = 0, g = 0, b = 0;
                switch (view) {
                case View::Dye: {
                    float d = clampf(m_dye[i], 0.f, 1.f);
                    float s = clampf(std::sqrt(m_u[i] * m_u[i] + m_v[i] * m_v[i]) * invU, 0.f, 1.f);
                    r = d * (0.25f + 0.75f * s);
                    g = d * (0.55f + 0.35f * s);
                    b = d * (0.95f);
                    b += (1.f - d) * 0.04f;
                    break;
                }
                case View::Speed: {
                    float s = clampf(std::sqrt(m_u[i] * m_u[i] + m_v[i] * m_v[i]) * invU, 0.f, 1.f);
                    r = clampf(s * 2.2f, 0.f, 1.f);
                    g = clampf(s * 1.4f - 0.25f, 0.f, 1.f);
                    b = clampf(0.35f * s + (s > 0.85f ? (s - 0.85f) * 4.f : 0.f), 0.f, 1.f);
                    break;
                }
                case View::Vorticity: {
                    float w = 0.f;
                    if (x > 0 && x < nx - 1 && y > 0 && y < ny - 1) {
                        w = (m_v[idx(x + 1, y)] - m_v[idx(x - 1, y)])
                          - (m_u[idx(x, y + 1)] - m_u[idx(x, y - 1)]);
                        w *= 18.f / std::max(params.u0, 1e-4f) * 0.05f;
                    }
                    float t = clampf(w, -1.f, 1.f);
                    if (t >= 0) { r = t; g = 0.15f * t; b = 0.06f; }
                    else        { b = -t; g = 0.25f * -t; r = 0.05f; }
                    r += 0.03f; g += 0.03f; b += 0.05f;
                    break;
                }
                case View::Density: {
                    float d = clampf((m_rho[i] - 0.97f) / 0.06f, 0.f, 1.f);
                    r = d; g = 0.3f + 0.5f * d; b = 1.f - d * 0.7f;
                    r *= 0.8f; g *= 0.8f; b *= 0.8f;
                    break;
                }
                }
                out[x] = qRgb(int(clampf(r, 0, 1) * 255),
                              int(clampf(g, 0, 1) * 255),
                              int(clampf(b, 0, 1) * 255));
            }
        });
    }

    int m_nx = 600, m_ny = 300;
    std::vector<float> m_f[Q], m_fp[Q];
    std::vector<float> m_u, m_v, m_rho, m_dye, m_dyeTmp;
    std::vector<uint8_t> m_solid;
    QImage m_image;
    QRect m_view;
    QTimer m_timer;
    QElapsedTimer m_statClock;
    float m_msPerStep = 0.f;
};

// -------------------------------------------------------------------- window

class MainWindow : public QWidget
{
    Q_OBJECT
public:
    MainWindow()
    {
        setWindowTitle("VortexStreet (D2Q9 lattice Boltzmann)");
        auto *canvas = new LbmCanvas;

        auto *grid = new QComboBox;
        grid->addItems({ "400 x 200", "600 x 300", "900 x 450" });
        grid->setCurrentIndex(1);

        auto *collision = new QComboBox;
        collision->addItems({ "BGK", "MRT" });
        collision->setCurrentIndex(1);

        auto *omega = new QDoubleSpinBox;
        omega->setRange(0.6, 1.99); omega->setDecimals(3); omega->setSingleStep(0.01);
        omega->setValue(canvas->params.omega);

        auto *u0 = new QDoubleSpinBox;
        u0->setRange(0.01, 0.14); u0->setDecimals(3); u0->setSingleStep(0.005);
        u0->setValue(canvas->params.u0);

        auto *steps = new QSpinBox;
        steps->setRange(1, 60);
        steps->setValue(canvas->params.stepsPerFrame);

        auto *view = new QComboBox;
        view->addItems({ "Dye", "Speed", "Vorticity", "Density" });

        auto *inject = new QCheckBox("Inject dye");
        inject->setChecked(true);

        auto *clearObs = new QPushButton("Clear obstacles");
        auto *reset = new QPushButton("Reset flow");
        auto *resetAll = new QPushButton("Reset flow + cylinder");
        auto *pause = new QPushButton("Pause");
        pause->setCheckable(true);

        auto *hint = new QLabel("LMB: draw obstacle\nRMB: erase");
        hint->setStyleSheet("color:#889");
        auto *stats = new QLabel;
        stats->setStyleSheet("color:#8fa;font-family:monospace");

        auto *form = new QFormLayout;
        form->addRow("Grid", grid);
        form->addRow("Collision", collision);
        form->addRow("Omega", omega);
        form->addRow("Inflow u0", u0);
        form->addRow("Steps/frame", steps);
        form->addRow("View", view);
        form->addRow(inject);
        form->addRow(clearObs);
        form->addRow(reset);
        form->addRow(resetAll);
        form->addRow(pause);
        form->addRow(hint);
        form->addRow(stats);

        auto *group = new QGroupBox("Parameters");
        group->setLayout(form);
        group->setFixedWidth(280);

        auto *layout = new QHBoxLayout(this);
        layout->addWidget(group);
        layout->addWidget(canvas, 1);

        connect(grid, &QComboBox::currentIndexChanged, this, [canvas](int i) {
            canvas->setGridPreset(i);
        });
        connect(collision, &QComboBox::currentIndexChanged, this, [canvas](int i) {
            canvas->params.collision = Collision(i);
        });
        connect(omega, &QDoubleSpinBox::valueChanged, this, [canvas](double v) { canvas->params.omega = float(v); });
        connect(u0, &QDoubleSpinBox::valueChanged, this, [canvas](double v) { canvas->params.u0 = float(v); });
        connect(steps, &QSpinBox::valueChanged, this, [canvas](int v) { canvas->params.stepsPerFrame = v; });
        connect(view, &QComboBox::currentIndexChanged, this, [canvas](int i) { canvas->params.view = View(i); });
        connect(inject, &QCheckBox::toggled, this, [canvas](bool b) { canvas->params.injectDye = b; });
        connect(clearObs, &QPushButton::clicked, this, [canvas] { canvas->clearObstacles(); });
        connect(reset, &QPushButton::clicked, this, [canvas] { canvas->reset(false); });
        connect(resetAll, &QPushButton::clicked, this, [canvas] { canvas->reset(true); });
        connect(pause, &QPushButton::toggled, this, [canvas, pause](bool b) {
            canvas->params.paused = b;
            pause->setText(b ? "Resume" : "Pause");
        });
        connect(canvas, &LbmCanvas::statsChanged, this, [canvas, stats] {
            stats->setText(QString("%1 ms/step\nRe ~ %2")
                           .arg(canvas->msPerStep(), 0, 'f', 2)
                           .arg(canvas->reynolds(), 0, 'f', 0));
        });

        resize(1260, 720);
    }
};

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    app.setStyle(QStyleFactory::create("Fusion"));
    QPalette pal;
    pal.setColor(QPalette::Window, QColor(37, 37, 42));
    pal.setColor(QPalette::WindowText, QColor(220, 220, 224));
    pal.setColor(QPalette::Base, QColor(28, 28, 32));
    pal.setColor(QPalette::Text, QColor(220, 220, 224));
    pal.setColor(QPalette::Button, QColor(48, 48, 54));
    pal.setColor(QPalette::ButtonText, QColor(220, 220, 224));
    pal.setColor(QPalette::Highlight, QColor(70, 120, 200));
    app.setPalette(pal);

    MainWindow w;
    w.show();
    return app.exec();
}

#include "main.moc"
