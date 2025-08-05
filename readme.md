# GeoCache: Kepler

## Problem
Ein Satelit, der sich auf einer kreisförmigen
Bahn um einen Planeten bewegt,
wurde von einem Gesteinsbrocken getroffen
und ist deshalb auf den Planeten gestürzt.

Der Satelit ist auf einer Bahn mit Radius
$R_{\text{Satelit}} =$ unterwegs.
Zum Zeitpunkt des Aufpralls hat sich der Satelit 
bei genau
TODO Grad Osten
über dem Äquator befunden und seine Geschwindigkeit hat
(von Norden aus gemessen), genau 45° in Richtung Osten gezeigt.

Durch den Aufprall wurde der Satelit abgebremst,
aber nicht abgelenkt, und seine Geschwindigkeit betrug
nach dem Aufprall noch
TODO%
seiner Anfangsgeschwindigkeit.

Der kreisförmige Planet dreht sich von Westen nach Osten,
hat einen Radius von
$R_{\text{Planet}} = $,
eine Tagesdauer von $24$ h,
und eine Mittlere Dichte von
$\rho_{\text{Planet}} = $.

Wo ist der Satelit auf den Planeten gestürzt?

## Solving the equations of motion

The Lagrangian of the system is given by
$$
    L
    :=
    \frac{m}{2}
    \left(
        \dot r^2 + r^2 \dot \theta^2
    \right)
    - V(r),
$$

where the Gravitational Potential is

$$
    V(r) = - \frac{k}{r},
$$

with

$$
    k := GMm.
$$

The Euler-Lagrange Equations are therefore given by

$$
    \begin{align}
        \frac{\mathrm d}{\mathrm d t}
        \left(
            m r^2 \dot \theta
        \right)
        &= 0, \\
        m \ddot r - mr \dot \theta^2 - \frac{\partial V}{\partial r} &= 0.
    \end{align}
$$

From the former, we may directly conclude
the _angular momentum_ $\ell$ to be conserved,
that is

$$
    \ell := m r^2 \dot \theta \equiv \mathrm{const}.
$$

The latter implies the conservation of total energy,
that is

$$
    E
    := \frac{1}{2} m \dot r^2
    + \frac{\ell^2}{2mr^2}
    + V(r) \equiv \mathrm{const},
$$

and hence, the Euler-Lagrange equations reduce to the
first order system

$$
\begin{equation}
    \begin{aligned}
        \dot \theta 
        &=
        \frac{\ell}{m r^2},\\
        \dot r
        &=
        \sqrt{
            \frac{2}{m}
            \left(
                E - V(r) - \frac{\ell^2}{2mr^2}
            \right)
        }
    \end{aligned}
\end{equation}
$$

Defining the _effective potential_ as

$$
    V_{\text{eff}}
    :=
    V(r) + \frac{\ell^2}{2mr^2},
$$

the first order system reduces to

$$
\begin{equation}
    \begin{aligned}
        \dot \theta 
        &=
        \frac{\ell}{m r^2},\\
        \dot r
        &=
        \sqrt{
            \frac{2}{m}
            \left(
                E - V_{\text{eff}}
            \right)
        }.
    \end{aligned}
\end{equation}
$$

