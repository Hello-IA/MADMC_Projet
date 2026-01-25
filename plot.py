import matplotlib.pyplot as plt

# Lorenz-efficient set (same for all omegas)
lorenz_points = [
    (4851, 11562),
    (5136, 11432),
    (5296, 11286),
]

lorenz_per_omega = {
    (1.01, 1): [
        (4851, 11562),  # P1
        (5136, 11432),
        (5296, 11286),
    ],
    (1.5, 1): [
        (5136, 11432),  # P1
        (4851, 11562),
        (5296, 11286),
    ],
    (2, 1): [
        (5296, 11286),  # P1
        (5136, 11432),
        (4851, 11562),
    ],
    (5, 1): [
        (5296, 11286),  # P1
        (5136, 11432),
        (4851, 11562),
    ],
    (10, 1): [
        (5296, 11286),  # P1
        (5136, 11432),
        (4851, 11562),
    ],
}

lorenz_frontier = [
    (4851, 11562),
    (5136, 11432),
    (5296, 11286),
]

# P1-selected Lorenz vector for each omega
omega_to_p1 = {
    (1.01, 1): (4851, 11562),
    (1.5, 1):  (5136, 11432),
    (2, 1):    (5296, 11286),
    (5, 1):    (5296, 11286),
    (10, 1):   (5296, 11286),
}


L1 = [p[0] for p in lorenz_points]
L2 = [p[1] for p in lorenz_points]
import matplotlib.pyplot as plt

plt.figure(figsize=(7,5))

# Plot Lorenz-efficient frontier
L1, L2 = zip(*lorenz_frontier)
plt.plot(L1, L2, "--o", color="black", label="Lorenz-efficient frontier")

# Colors for each omega
colors = ["red", "blue", "green", "purple", "orange"]

for (omega, color) in zip(lorenz_per_omega.keys(), colors):
    p1_point = lorenz_per_omega[omega][0]
    plt.scatter(*p1_point, s=120, facecolors="none",
                edgecolors=color, linewidths=2,
                label=f"P1 for ω={omega}")

plt.xlabel(r"$L_1$ (minimum objective)")
plt.ylabel(r"$L_2$ (sum of objectives)")
plt.title("Lorenz-efficient frontier and P1 selections for different ω")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
# plt.show()












import matplotlib.pyplot as plt

# Lorenz vectors per omega in the order you printed (Solution 0 = P1)
lorenz_per_omega = {
    (1.01, 1): [(4851, 11562), (5136, 11432), (5296, 11286)],
    (1.5, 1):  [(5136, 11432), (4851, 11562), (5296, 11286)],
    (2, 1):    [(5296, 11286), (5136, 11432), (4851, 11562)],
    (5, 1):    [(5296, 11286), (5136, 11432), (4851, 11562)],
    (10, 1):   [(5296, 11286), (5136, 11432), (4851, 11562)],
}

# Optional: your OWA values (same order as above). Used only for labels if you want.
owa_per_omega = {
    (1.01, 1): [11610.51, 11483.36, 11338.96],
    (1.5, 1):  [14000.00, 13987.50, 13934.00],
    (2, 1):    [16582.00, 16568.00, 16413.00],
    (5, 1):    [32470.00, 31976.00, 30966.00],
    (10, 1):   [58950.00, 57656.00, 55221.00],
}

# ---- Plot ----
plt.figure(figsize=(8, 6))

# Plot the Lorenz frontier once (order-independent)
frontier = sorted({pt for seq in lorenz_per_omega.values() for pt in seq})
fx, fy = zip(*frontier)
plt.plot(fx, fy, "--", color="gray", alpha=0.6, label="Lorenz-efficient frontier")

# Colors for each omega
omegas = list(lorenz_per_omega.keys())
colors = ["red", "blue", "green", "purple", "orange"]

for omega, color in zip(omegas, colors):
    pts = lorenz_per_omega[omega]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]

    # PL points (all enumerated Lorenz vectors for this omega)
    plt.scatter(xs, ys, s=70, color=color, alpha=0.35)

    # Connect in enumeration order (optional but useful)
    plt.plot(xs, ys, "-", color=color, alpha=0.35)

    # Highlight P1 (first point)
    p1x, p1y = pts[0]
    plt.scatter(p1x, p1y, s=170, facecolors="none", edgecolors=color,
                linewidths=2.5, label=f"P1 for ω={omega}")

    # Optional label at P1
    plt.text(p1x + 12, p1y - 25, f"ω={omega}", fontsize=9)

# Axis labels and title
plt.xlabel(r"$L_1$ (minimum objective)")
plt.ylabel(r"$L_2$ (sum of objectives)")
plt.title("P1 selections vs PL enumerations in Lorenz space")
plt.grid(alpha=0.25)
plt.legend()
plt.tight_layout()
plt.show()


