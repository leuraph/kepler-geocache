# GeoCache: Kepler

## Problem


Ein Satelit, der sich auf einer kreisförmigen
Bahn um einen Planeten bewegt,
ist mit einem unbekannten Objekt kollidiert
und deshalb auf den Planeten gestürzt.

Der kugelförmige Planet dreht sich
mit einer Tagesdauer von $24~\text{h}$
von Westen nach Osten,
hat einen Radius von $R_{\text{P}} = 5700~\text{km}$,
und eine Dichte von
$\rho_{\text{P}} = 3.344~\frac{\text{g}}{\text{cm}^3}$.

Vor der Kollision war der Satelit
auf einer kreisförmigen Bahn mit Radius
$R_{\text{S}} =2 R_{\text{P}}$ unterwegs.
Zum Zeitpunkt der Kollision mit dem unbekannten Objekt
hat sich der Satelit bei genau
$-5.287518$
Grad Osten über dem Äquator befunden und seine Geschwindigkeit hat
(von Norden aus gemessen), genau
$\frac{\pi}{8}$ (in Bogenmass)
in Richtung Osten gezeigt.

Durch den Aufprall wurde der Satelit **abgebremst**,
aber **nicht abgelenkt**, und seine Geschwindigkeit betrug
nach dem Aufprall noch
$53.283236\%$
von seiner Anfangsgeschwindigkeit.

Wo ist der Satelit auf den Planeten gestürzt?

## Computing Initial Conditions

Zuerst bestimmen wir die Masse des Planeten mit
```math
M_{\text P} = \rho_{\text P} \frac{4\pi}{3} R_{\text{P}}^3.
```

Um die Bahngeschwindigkeit $v_{\text{S}}$ des Sateliten
vor der Kollision zu bestimmen, setzen wir die Zentripetalkraft
mit der Gravitationskraft gleich, i.e.

```math
\frac{Gm_{\text{S}}M_{\text{P}}}{R_{\text{S}}^2}
=
\frac{m_{\text{S}}v_{\text{S}}^2}{R_{\text{S}}},
```

wobei sich die Masse des Sateliten $m_{\text S}$ sowieso rausstreicht.
Wir finden also

```math
v_{\text{S}}
= \sqrt{\frac{GM_{\text{P}}}{R_{\text{S}}}}
= \sqrt{\frac{GM_{\text{P}}}{2 R_{\text{P}}}}.
```

Wenn $\gamma = 0.53283236$,
dann ist die Geschwindigkeit des Sateliten nach der Kollision
gegenen durch $\gamma v_{\text S}$.
Der Geschwindigkeitsvektor unmittelbar nach der Kollision ist demnach

```math
v = 
(0, \, 
\sin(\alpha) \gamma v_\text{S}, \,
\cos(\alpha) \gamma v_\text{S}),
```
wobei $\alpha = \frac{\pi}{8}$.
Die Position des Sateliten unmittelbar nach der Kollision ist
```math
u = (R_{\text S}, 0, 0).
```

## Equations of motion (Newtonian Mechanics)


We assume the spherical planet to be centered at the origin.
Then, if $u = (x,y,z) \in \mathbb{R}^3$
denotes the satelite's position in standard euclidean coordinates,
the Gravitational force acting on it is given by
$F(u) = - \frac{GMm}{\|u\|^3} u$.
The correspondind Equation of motion (EOM)
is then given by
```math
\ddot u = - \frac{GM}{\|u\|^3} u.
```
Writing
$k:= GM$, $r := \sqrt{x^2 + y^2 + z^2}$, and
$Y = (x,y,z,\dot x, \dot y, \dot z)$
yields the ODE
```math
\dot Y
=
(\dot x, \dot y, \dot z, \ddot x, \ddot y, \ddot z)
=
(Y^3, Y^4, Y^5,
-\frac{k}{r^3} Y^0,
-\frac{k}{r^3} Y^1,
-\frac{k}{r^3} Y^2).
```
Given the initial conditions from above, i.e.

```math
Y_0 =
(R_{\text{S}}, 0, 0, 0,
\sin(\alpha) \gamma v_\text{S},
\cos(\alpha) \gamma v_\text{S}),
```

where $\alpha = \frac{\pi}{8}$
and $\gamma = 0.53283236$,
the ODE may be numerically integrated over
$(0, t_{\text{end}})$, $t_{\text{end}} > 0$.
To find the coordinates of impact,
one may define a function
$F : (0, \infty) \to \mathbb R$,

```math
F(t_{\text{end}})
:=
R(\tilde u (t_{\text{end}})) - R_{\text{P}},
```

where $\tilde u$ is the numerically integrated solution of the ODE,
$R(u) := \sqrt{x^2 + y^2 + z^2}$,
and $R_{\text P}$ denotes the radius of the planet,
and then apply a root finding algorithm to $F$.
Then, the coordinates of impact may first be transformed
to spherical coordinates
$(r_{\text{hit}}, \theta_{\text{hit}}, \phi_{\text{hit}})$,
where $r$ denotes the radial distance to the origin,
$\theta$ is the polar angle measured from the positive $z$-axis,
and $\phi$ the azimuthal angle measured from the positive $x$-axis
(see [wiki](https://en.wikipedia.org/wiki/Spherical_coordinate_systemhttps://en.wikipedia.org/wiki/Spherical_coordinate_system)).
Finally, the coordinates of impact may be calculated with

```math
(N, E)
=
\bigg(
    \frac{360}{2\pi}\bigg(\frac{\pi}{2} - \theta_{\text{hit}}\bigg),\,
    E_{\text{initial}} + \frac{360}{2\pi}\bigg( \phi_{\text{hit}} - t_{\text{hit}} \omega_{\text{P}} \bigg)
\bigg),
```

where
$E_{\text{initial}} = −5.287518$,
and $\omega_{\text P}$ denotes the planet's angular velocity
($\omega = 2\pi f = \frac{2\pi}{T}$).

## Equations of motion (Kepler)


In the following, we drop the indices $\text S$ and $\text P$.
Due to the conservation of angular momentum,
the motion of the satellite takes place in a plane.
For convenience, we work in polar coordinates $(r, \theta)$.
The Lagrangian of the system is given by

```math
    L
    :=
    \frac{m}{2}
    \left(
        \dot r^2 + r^2 \dot \theta^2
    \right)
    - V(r),
```

where the Gravitational Potential $V$ is given by

```math
    V(r) = - \frac{GMm}{r},
```

The Euler-Lagrange Equations are given by

```math
    \begin{align}
        \frac{\mathrm d}{\mathrm d t}
        \left(
            m r^2 \dot \theta
        \right)
        &= 0, \\
        m \ddot r - mr \dot \theta^2 - \frac{\partial V}{\partial r} &= 0.
    \end{align}
```

From the former, we may directly conclude
the _angular momentum_ $\ell$ to be conserved,
that is

```math
    \ell := m r^2 \dot \theta \equiv \mathrm{const}.
```

It can be shown that the total energy $E$
of the system is conserved
(either using the EOM of the coordinate $r$
or by using $\mathrm d L/ \mathrm d t = 0$),
that is
```math
    E
    := \frac{1}{2} m \dot r^2
    + \frac{\ell^2}{2mr^2}
    + V(r) \equiv \mathrm{const}.
```

Hence, the Euler-Lagrange equations reduce to the
first order system

```math
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
```

Defining the _effective potential_ as

```math
    V_{\text{eff}}
    :=
    V(r) + \frac{\ell^2}{2mr^2},
```

the first order system reduces to

```math
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
```

