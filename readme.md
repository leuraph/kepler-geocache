# GeoCache: Kepler

## Problem
Ein Satelit, der sich auf einer kreisförmigen
Bahn um einen Planeten bewegt,
wurde von einem Gesteinsbrocken getroffen
und ist deshalb auf den Planeten gestürzt.

Der Satelit ist auf einer Bahn mit Radius
$R_{\text{S}} =2 R_{\text{P}}$ unterwegs.
Zum Zeitpunkt der Kollision mit dem unbekannten Objekt
hat sich der Satelit bei genau
$0$
Grad Osten über dem Äquator befunden und seine Geschwindigkeit hat
(von Norden aus gemessen), genau 45° in Richtung Osten gezeigt.

Durch den Aufprall wurde der Satelit abgebremst,
aber nicht abgelenkt, und seine Geschwindigkeit betrug
nach dem Aufprall noch
$\alpha = 0.5$
seiner Anfangsgeschwindigkeit.

Der kreisförmige Planet dreht sich von Westen nach Osten,
hat einen Radius von
$R_{\text{P}} = 5700$ km,
eine Tagesdauer von $24~\text{h}$,
und eine Mittlere Dichte von
$\rho_{\text{P}} = 3.344~\frac{\text{g}}{\text{cm}^3}$.

Wo ist der Satelit auf den Planeten gestürzt?

## Computing Initial Conditions

Um die Bahngeschwindigkeit $v_{\text{S}}$ des Sateliten
vor der Kollision zu bestimmen, setzen wir die Zentripetalkraft
mit der Gravitationskraft gleich, i.e.
$$
\frac{Gm_{\text{S}}M_{\text{P}}}{R_{\text{S}}^2}
=
\frac{m_{\text{S}}v_{\text{S}}^2}{R_{\text{S}}}.
$$
Wir finden also
$$
v_{\text{S}} = \sqrt{\frac{GM_{\text{P}}}{R_{\text{S}}}}.
$$

## Equations of motion

In the following, we drop the indices $\text S$ and $\text P$.
Due to the conservation of angular momentum,
the motion of the satellite takes place in a plane.
For convenience, we work in polar coordinates $(r, \theta)$.
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

where the Gravitational Potential $V$ is given by

$$
    V(r) = - \frac{GMm}{r},
$$

The Euler-Lagrange Equations are given by

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

It can be shown that the total energy $E$
of the system is conserved
(either using the EOM of the coordinate $r$
or by using $\mathrm d L/ \mathrm d t = 0$),
that is
$$
    E
    := \frac{1}{2} m \dot r^2
    + \frac{\ell^2}{2mr^2}
    + V(r) \equiv \mathrm{const}.
$$

Hence, the Euler-Lagrange equations reduce to the
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

