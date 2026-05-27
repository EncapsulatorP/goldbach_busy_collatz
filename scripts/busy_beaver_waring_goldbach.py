#!/usr/bin/env python3
"""
Busy Beaver space-time generator with a Waring-Goldbach-style tableau encoding.

This is not an exact Busy Beaver solver for arbitrary state counts. Instead, it
simulates small built-in Busy Beaver machines and lifts their computation
history into a prime-power tableau with the following premises:

1. Each visible spacetime cell (t, x) receives a unique prime p_{t,x}.
2. p_{t,x} mod 30 encodes the local symbol/head class:
       1  -> blank, no head
       7  -> one,   no head
       11 -> blank, head present
       13 -> one,   head present
3. The exponent k_{t,x} encodes time-scale depth, head/state prominence, and
   the outgoing transition type from that head position.
4. Each time row becomes a Waring-Goldbach-like mass:
       N_t = sum_x p_{t,x}^{k_{t,x}}
5. Each executed machine step emits a 4-prime "Goldbach plaquette" over
   (t, x-1), (t, x), (t, x+1), and the next head position at time t+1.

The result is a computable toy encoding of Busy Beaver spacetime, not a proof
mechanism for the Busy Beaver function or for Goldbach/Waring-Goldbach claims.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import pandas as pd
from sympy import nextprime

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


MODULUS = 30
RESIDUE_MAP = {
    (0, False): 1,
    (1, False): 7,
    (0, True): 11,
    (1, True): 13,
}
TRANSITION_RANK = {
    "L": 1,
    "R": 2,
    "HALT": 3,
    "STOP": 4,
}

# Shell labels: residue mod 30 → (description, plot colour)
# Each residue class is one "prime-boost shell" in the space-time tableau.
SHELL_META: dict[int, tuple[str, str]] = {
    1:  ("blank  / no head",    "#4a9eff"),   # blue  — background tape
    7:  ("symbol=1 / no head",  "#ff7c43"),   # amber — written cells
    11: ("blank  + head",       "#a8e6a3"),   # green — head on blank
    13: ("symbol=1 + head",     "#ffd700"),   # gold  — head on symbol=1
}


@dataclass(frozen=True)
class Transition:
    write: int
    move: int
    next_state: str


@dataclass(frozen=True)
class MachineSpec:
    name: str
    states: tuple[str, ...]
    start_state: str
    halt_state: str
    transitions: dict[tuple[str, int], Transition]


@dataclass(frozen=True)
class Configuration:
    time: int
    state: str
    head: int
    tape: dict[int, int]
    halted: bool


def builtin_machines() -> dict[str, MachineSpec]:
    """
    Small verified 2-symbol Busy Beaver machines.

    The built-in set is intentionally conservative: these are compact machines
    that produce interesting halting space-time histories quickly.
    """
    return {
        "bb2": MachineSpec(
            name="bb2",
            states=("A", "B"),
            start_state="A",
            halt_state="H",
            transitions={
                ("A", 0): Transition(write=1, move=1, next_state="B"),
                ("A", 1): Transition(write=1, move=-1, next_state="B"),
                ("B", 0): Transition(write=1, move=-1, next_state="A"),
                ("B", 1): Transition(write=1, move=1, next_state="H"),
            },
        ),
        "bb3": MachineSpec(
            name="bb3",
            states=("A", "B", "C"),
            start_state="A",
            halt_state="H",
            transitions={
                ("A", 0): Transition(write=1, move=1, next_state="B"),
                ("A", 1): Transition(write=1, move=1, next_state="H"),
                ("B", 0): Transition(write=0, move=1, next_state="C"),
                ("B", 1): Transition(write=1, move=1, next_state="B"),
                ("C", 0): Transition(write=1, move=-1, next_state="C"),
                ("C", 1): Transition(write=1, move=-1, next_state="A"),
            },
        ),
        # BB(4) champion: 107 steps, 13 ones (Rado 1962 / Brady 1983)
        # This is the verified 4-state 2-symbol busy beaver.
        "bb4": MachineSpec(
            name="bb4",
            states=("A", "B", "C", "D"),
            start_state="A",
            halt_state="H",
            transitions={
                ("A", 0): Transition(write=1, move=1,  next_state="B"),
                ("A", 1): Transition(write=1, move=-1, next_state="B"),
                ("B", 0): Transition(write=1, move=-1, next_state="A"),
                ("B", 1): Transition(write=0, move=-1, next_state="C"),
                ("C", 0): Transition(write=1, move=1,  next_state="H"),
                ("C", 1): Transition(write=1, move=-1, next_state="D"),
                ("D", 0): Transition(write=1, move=1,  next_state="D"),
                ("D", 1): Transition(write=0, move=1,  next_state="A"),
            },
        ),
    }


def move_label(move: int) -> str:
    if move < 0:
        return "L"
    if move > 0:
        return "R"
    return "S"


def transition_code(transition: Transition) -> str:
    if transition.next_state == "H":
        return "HALT"
    return move_label(transition.move)


def ensure_parent_dir(path: str) -> None:
    parent = Path(path).expanduser().resolve().parent
    parent.mkdir(parents=True, exist_ok=True)


def sparse_write(tape: dict[int, int], position: int, value: int) -> None:
    if value:
        tape[position] = value
    elif position in tape:
        del tape[position]


def simulate_machine(machine: MachineSpec, max_steps: int) -> list[Configuration]:
    tape: dict[int, int] = {}
    head = 0
    state = machine.start_state
    history = [
        Configuration(
            time=0,
            state=state,
            head=head,
            tape=dict(tape),
            halted=False,
        )
    ]

    for step in range(1, max_steps + 1):
        read = int(tape.get(head, 0))
        transition = machine.transitions.get((state, read))
        if transition is None:
            raise ValueError(f"Missing transition for state={state}, symbol={read}")

        sparse_write(tape, head, transition.write)
        head += transition.move
        state = transition.next_state
        halted = state == machine.halt_state

        history.append(
            Configuration(
                time=step,
                state=state,
                head=head,
                tape=dict(tape),
                halted=halted,
            )
        )

        if halted:
            break

    return history


def history_span(history: list[Configuration], padding: int) -> list[int]:
    positions = {cfg.head for cfg in history}
    for cfg in history:
        positions.update(cfg.tape.keys())

    lo = min(positions) - padding
    hi = max(positions) + padding
    return list(range(lo, hi + 1))


def scale_level(time_index: int) -> int:
    if time_index <= 0:
        return 1
    return 1 + int(math.log2(time_index + 1))


def outgoing_transition(machine: MachineSpec, cfg: Configuration) -> Transition | None:
    if cfg.halted or cfg.state == machine.halt_state:
        return None
    read = int(cfg.tape.get(cfg.head, 0))
    return machine.transitions[(cfg.state, read)]


def next_prime_with_residue(start_after: int, residue: int, modulus: int) -> int:
    candidate = max(1, start_after)
    while True:
        candidate = int(nextprime(candidate))
        if candidate % modulus == residue:
            return candidate


def state_rank_map(machine: MachineSpec) -> dict[str, int]:
    ranks = {state: idx + 1 for idx, state in enumerate(machine.states)}
    ranks[machine.halt_state] = len(machine.states) + 1
    return ranks


def build_cell_records(
    history: list[Configuration],
    machine: MachineSpec,
    span: list[int],
    modulus: int = MODULUS,
) -> tuple[list[dict], dict[tuple[int, int], dict]]:
    residue_cursor = {residue: 1 for residue in set(RESIDUE_MAP.values())}
    state_ranks = state_rank_map(machine)
    records: list[dict] = []
    lookup: dict[tuple[int, int], dict] = {}

    for cfg in history:
        transition = outgoing_transition(machine, cfg)
        t_scale = scale_level(cfg.time)

        for x in span:
            symbol = int(cfg.tape.get(x, 0))
            is_head = x == cfg.head
            residue = RESIDUE_MAP[(symbol, is_head)]
            prime = next_prime_with_residue(residue_cursor[residue], residue, modulus)
            residue_cursor[residue] = prime

            state_here = cfg.state if is_head else ""
            transition_here = transition_code(transition) if is_head and transition else ""
            state_rank = state_ranks.get(cfg.state, 0) if is_head else 0
            transition_rank = TRANSITION_RANK.get(transition_here, 0)
            exponent = 2 + t_scale + symbol + int(is_head) + state_rank + transition_rank
            term_value = prime**exponent
            term_log10 = math.log10(term_value) if term_value > 0 else 0.0

            record = {
                "time": cfg.time,
                "x": x,
                "symbol": symbol,
                "is_head": int(is_head),
                "state": state_here,
                "transition_code": transition_here,
                "residue": residue,
                "prime": prime,
                "exponent": exponent,
                "scale_level": t_scale,
                "term_value": term_value,
                "term_log10": term_log10,
            }
            records.append(record)
            lookup[(cfg.time, x)] = record

    return records, lookup


def build_row_records(
    history: list[Configuration],
    span: list[int],
    cell_lookup: dict[tuple[int, int], dict],
) -> list[dict]:
    rows: list[dict] = []

    for cfg in history:
        row_cells = [cell_lookup[(cfg.time, x)] for x in span]
        row_mass = sum(cell["term_value"] for cell in row_cells)
        ones_count = sum(int(cfg.tape.get(x, 0)) for x in span)
        active_positions = sorted(set(cfg.tape.keys()) | {cfg.head})
        active_left = active_positions[0]
        active_right = active_positions[-1]
        support_width = active_right - active_left + 1

        rows.append(
            {
                "time": cfg.time,
                "state": cfg.state,
                "head": cfg.head,
                "halted": int(cfg.halted),
                "ones_count": ones_count,
                "active_left": active_left,
                "active_right": active_right,
                "support_width": support_width,
                "row_mass": row_mass,
                "row_mass_log10": math.log10(row_mass) if row_mass > 0 else 0.0,
            }
        )

    return rows


def build_transition_plaquettes(
    history: list[Configuration],
    machine: MachineSpec,
    cell_lookup: dict[tuple[int, int], dict],
    modulus: int = MODULUS,
) -> list[dict]:
    state_ranks = state_rank_map(machine)
    tiles: list[dict] = []

    for idx in range(len(history) - 1):
        cfg = history[idx]
        next_cfg = history[idx + 1]
        transition = outgoing_transition(machine, cfg)
        if transition is None:
            continue

        x0 = cfg.head
        x1 = next_cfg.head
        left = cell_lookup[(cfg.time, x0 - 1)]
        center = cell_lookup[(cfg.time, x0)]
        right = cell_lookup[(cfg.time, x0 + 1)]
        next_head = cell_lookup[(next_cfg.time, x1)]

        plaquette_sum = left["prime"] + center["prime"] + right["prime"] + next_head["prime"]
        plaquette_power_sum = (
            left["prime"]
            + center["prime"] ** 2
            + right["prime"] ** 3
            + next_head["prime"] ** 4
        )
        rule_code = transition_code(transition)
        rule_residue = (
            left["residue"]
            + 2 * center["residue"]
            + 3 * right["residue"]
            + 5 * next_head["residue"]
            + 7 * state_ranks[cfg.state]
            + 11 * transition.write
            + 13 * TRANSITION_RANK[rule_code]
        ) % modulus

        tiles.append(
            {
                "time": cfg.time,
                "source_state": cfg.state,
                "read_symbol": int(cfg.tape.get(cfg.head, 0)),
                "write_symbol": transition.write,
                "move": move_label(transition.move),
                "next_state": transition.next_state,
                "source_head": x0,
                "next_head": x1,
                "local_word": f"{left['symbol']}{center['symbol']}{right['symbol']}",
                "plaquette_residues": (
                    f"{left['residue']},{center['residue']},"
                    f"{right['residue']},{next_head['residue']}"
                ),
                "plaquette_sum": plaquette_sum,
                "plaquette_sum_mod": plaquette_sum % modulus,
                "plaquette_power_sum": plaquette_power_sum,
                "plaquette_power_sum_mod": plaquette_power_sum % modulus,
                "rule_residue": rule_residue,
                "rule_code": f"{cfg.state}{int(cfg.tape.get(cfg.head, 0))}"
                f"->{transition.write}{move_label(transition.move)}{transition.next_state}",
            }
        )

    return tiles


def stringify_large_int_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column in out.columns:
            out[column] = out[column].map(str)
    return out


def plot_tableau(
    cell_df: pd.DataFrame,
    row_df: pd.DataFrame,
    machine: MachineSpec,
    out_path: str,
) -> None:
    ensure_parent_dir(out_path)

    pivot = (
        cell_df.pivot(index="time", columns="x", values="symbol")
        .sort_index(axis=0)
        .sort_index(axis=1)
    )
    head_rows = cell_df[cell_df["is_head"] == 1][["time", "x", "state"]].sort_values("time")

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(14, 11),
        gridspec_kw={"height_ratios": [3.0, 1.4, 1.2]},
        constrained_layout=True,
    )
    ax_grid, ax_mass, ax_head = axes

    cmap = ListedColormap(["#f4f1de", "#e07a5f"])
    im = ax_grid.imshow(pivot.values, aspect="auto", interpolation="nearest", cmap=cmap)
    ax_grid.set_title(f"{machine.name} spacetime tableau")
    ax_grid.set_ylabel("time")
    ax_grid.set_xlabel("tape position")
    ax_grid.set_xticks(range(len(pivot.columns)))
    ax_grid.set_xticklabels([str(x) for x in pivot.columns], rotation=0)

    x_to_idx = {x: idx for idx, x in enumerate(pivot.columns)}
    head_x = [x_to_idx[int(x)] for x in head_rows["x"]]
    head_y = [int(t) for t in head_rows["time"]]
    ax_grid.plot(head_x, head_y, color="#1d3557", linewidth=1.8, marker="o", markersize=4)

    for row in head_rows.itertuples(index=False):
        ax_grid.text(
            x_to_idx[int(row.x)] + 0.08,
            int(row.time) + 0.12,
            row.state,
            color="#081c15",
            fontsize=8,
        )

    cbar = fig.colorbar(im, ax=ax_grid, fraction=0.03, pad=0.02)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["0", "1"])

    ax_mass.plot(row_df["time"], row_df["row_mass_log10"], color="#264653", linewidth=2.0)
    ax_mass.set_title("Row Waring-Goldbach mass")
    ax_mass.set_ylabel("log10(row mass)")
    ax_mass.set_xlabel("time")
    ax_mass.grid(alpha=0.25, linewidth=0.6)

    ax_head.plot(row_df["time"], row_df["head"], color="#2a9d8f", linewidth=2.0, label="head x")
    ax_head.plot(row_df["time"], row_df["ones_count"], color="#bc6c25", linewidth=2.0, label="ones")
    ax_head.set_title("Head position and ones count")
    ax_head.set_xlabel("time")
    ax_head.grid(alpha=0.25, linewidth=0.6)
    ax_head.legend(loc="best")

    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_shell_tableau(
    cell_df: pd.DataFrame,
    row_df: pd.DataFrame,
    machine: MachineSpec,
    out_path: str,
) -> None:
    """
    Four-panel shell-shaped space-time diagram.

    The encoding maps every (time, x) cell to a prime p ≡ r (mod 30) where r
    is one of four residue classes called "prime-boost shells":

        residue  1  →  blank tape cell, no head  (background)
        residue  7  →  symbol-1 cell, no head    (written mark)
        residue 11  →  blank cell WITH head       (reading blank)
        residue 13  →  symbol-1 cell WITH head    (reading one)

    The "prime boost" of a cell is its exponent k:  the cell contributes
    p^k to the row's Waring-Goldbach mass.  Larger k = deeper computational
    involvement (head present, deeper state, higher transition rank).

    Panel A — Classic spacetime tape with head path and exponent overlay
    Panel B — Shell-shaped concentric-ring diagram
               Each ring = one prime-boost shell (residue class).
               Cells are arranged by angle (tape x position, wrapping)
               and radially by time (inner = early, outer = late).
               Colour = exponent k (prime boost magnitude).
    Panel C — Prime-boost (exponent) evolution per shell
               X = time step, Y = exponent k, one line per shell.
               Shows how each shell's boost grows as the machine runs.
    Panel D — Waring-Goldbach row-mass spiral (log-polar)
               The row mass N_t = Σ_x p^k grows with time.
               Plotted in polar: angle ∝ time, radius = log10(row mass).
               This literally spirals outward like a nautilus shell as
               the machine's computation accumulates mass.
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.collections import LineCollection

    ensure_parent_dir(out_path)

    shells_in_data = sorted(cell_df["residue"].unique())
    max_time = int(cell_df["time"].max())
    max_exp  = int(cell_df["exponent"].max())
    exp_norm = Normalize(vmin=int(cell_df["exponent"].min()), vmax=max_exp)

    pivot = (
        cell_df.pivot(index="time", columns="x", values="symbol")
        .sort_index(axis=0).sort_index(axis=1)
    )
    xs_all = sorted(pivot.columns)
    x_min, x_max = xs_all[0], xs_all[-1]
    x_span = max(1, x_max - x_min)

    fig, axes = plt.subplots(2, 2, figsize=(18, 14), constrained_layout=True)
    fig.patch.set_facecolor("#0d0d14")
    for ax in axes.ravel():
        ax.set_facecolor("#0d0d14")
        for spine in ax.spines.values():
            spine.set_edgecolor("#333340")
        ax.tick_params(colors="#b0b0c0")
        ax.title.set_color("#e8e8f0")
        ax.xaxis.label.set_color("#b0b0c0")
        ax.yaxis.label.set_color("#b0b0c0")

    ax_tape, ax_rings, ax_boost, ax_spiral = axes.ravel()

    # ── Panel A: classic tape + head path + exponent overlay ────────────────
    tape_cmap = plt.get_cmap("Greys_r")
    im_tape = ax_tape.imshow(
        pivot.values, aspect="auto", interpolation="nearest",
        cmap=tape_cmap, vmin=0, vmax=1,
    )
    head_rows = (
        cell_df[cell_df["is_head"] == 1][["time", "x", "state", "exponent"]]
        .sort_values("time")
    )
    x_to_col = {x: i for i, x in enumerate(xs_all)}
    head_xs  = [x_to_col[int(r.x)] for r in head_rows.itertuples()]
    head_ts  = [int(r.time) for r in head_rows.itertuples()]
    head_ks  = [int(r.exponent) for r in head_rows.itertuples()]

    ax_tape.plot(head_xs, head_ts, color="#ffd700", lw=2, alpha=0.9, zorder=3)
    sc = ax_tape.scatter(head_xs, head_ts, c=head_ks, cmap="plasma",
                         norm=exp_norm, s=60, zorder=4, edgecolors="#000000", lw=0.4)
    for i, row in enumerate(head_rows.itertuples(index=False)):
        ax_tape.text(
            x_to_col[int(row.x)] + 0.1, int(row.time) + 0.15,
            row.state, color="#ffffff", fontsize=7, zorder=5,
        )
    cb_a = fig.colorbar(sc, ax=ax_tape, fraction=0.04, pad=0.03)
    cb_a.set_label("exponent k  (prime boost)", color="#b0b0c0")
    plt.setp(cb_a.ax.yaxis.get_ticklabels(), color="#b0b0c0")
    cb_a.ax.yaxis.set_tick_params(color="#b0b0c0")

    ax_tape.set_title(
        f"A — {machine.name.upper()} Spacetime Tableau\n"
        "White = symbol 1,  black = blank.  "
        "Gold path = head trajectory,  dot colour = exponent k (prime boost).",
        pad=6,
    )
    ax_tape.set_xlabel("Tape position x")
    ax_tape.set_ylabel("Time step t")
    ax_tape.set_xticks(range(len(xs_all)))
    ax_tape.set_xticklabels([str(x) for x in xs_all], rotation=0)

    # ── Panel B: shell-shaped concentric-ring diagram ────────────────────────
    ring_w = 1.6
    ring_g = 0.5
    cmap_exp = plt.get_cmap("plasma")
    rng = __import__("numpy").random.default_rng(3)

    ax_rings.set_aspect("equal")

    for shell_idx, residue in enumerate(sorted(SHELL_META)):
        label, base_color = SHELL_META[residue]
        sub = cell_df[cell_df["residue"] == residue].sort_values(["time", "x"])
        if sub.empty:
            continue

        r_in  = shell_idx * (ring_w + ring_g)
        r_out = r_in + ring_w

        # angle: tape position x maps to [0, 2π], wrapping naturally
        angles = (
            ((sub["x"].values - x_min) / x_span * 2 * math.pi)
            + sub["time"].values / max(1, max_time) * math.pi * 0.3  # slight time twist
        ) % (2 * math.pi)
        r_vals = rng.uniform(r_in + 0.06, r_out - 0.06, len(sub))

        xs_r = r_vals * __import__("numpy").cos(angles)
        ys_r = r_vals * __import__("numpy").sin(angles)

        colors_r = cmap_exp(exp_norm(sub["exponent"].values))
        ax_rings.scatter(xs_r, ys_r, c=colors_r, s=18, alpha=0.82,
                         linewidths=0, zorder=2)

        # Ring boundary arc
        theta_arc = __import__("numpy").linspace(0, 2 * math.pi, 360)
        ax_rings.plot(
            r_out * __import__("numpy").cos(theta_arc),
            r_out * __import__("numpy").sin(theta_arc),
            color=base_color, lw=1.2, alpha=0.5, zorder=3,
        )
        # Shell label at 3 o'clock
        ax_rings.text(
            r_out + 0.15, 0,
            f"r≡{residue} mod 30\n{label}",
            va="center", ha="left", color=base_color,
            fontsize=7, alpha=0.9,
        )

    sm_b = plt.cm.ScalarMappable(cmap=cmap_exp, norm=exp_norm)
    sm_b.set_array([])
    cb_b = fig.colorbar(sm_b, ax=ax_rings, fraction=0.04, pad=0.04)
    cb_b.set_label("exponent k  (prime boost)", color="#b0b0c0")
    plt.setp(cb_b.ax.yaxis.get_ticklabels(), color="#b0b0c0")
    cb_b.ax.yaxis.set_tick_params(color="#b0b0c0")

    ax_rings.set_title(
        "B — Prime-Boost Shell Cross-Section\n"
        "Concentric rings = residue shells (mod 30 class of assigned prime).\n"
        "Angle = tape-x position,  inner→outer = earlier→later time,  "
        "colour = exponent k.",
        pad=6,
    )
    ax_rings.axis("off")

    # ── Panel C: exponent (prime boost) evolution per shell ──────────────────
    time_steps = sorted(cell_df["time"].unique())
    for residue in sorted(SHELL_META):
        label, color = SHELL_META[residue]
        sub = cell_df[cell_df["residue"] == residue]
        mean_k = sub.groupby("time")["exponent"].mean().reindex(time_steps)
        max_k  = sub.groupby("time")["exponent"].max().reindex(time_steps)
        ax_boost.plot(time_steps, mean_k.values, color=color, lw=2,
                      label=f"r≡{residue}  {label}")
        ax_boost.fill_between(time_steps, mean_k.values, max_k.values,
                              color=color, alpha=0.15)

    ax_boost.set_title(
        "C — Prime-Boost (Exponent k) Evolution per Shell\n"
        "Solid line = mean exponent,  shaded band = mean → max.\n"
        "Spikes mark when the head enters a shell (exponent jumps).",
        pad=6,
    )
    ax_boost.set_xlabel(
        "Time step t\n"
        "Exponent k = 2 + scale_level + symbol + is_head + state_rank + transition_rank"
    )
    ax_boost.set_ylabel("Exponent k  (prime boost magnitude)")
    ax_boost.legend(fontsize=8, facecolor="#1a1a22", edgecolor="#444455",
                    labelcolor="#e0e0e8")
    ax_boost.grid(True, alpha=0.15, color="#444455")

    # ── Panel D: Waring-Goldbach row-mass log-polar spiral ───────────────────
    times   = row_df["time"].values
    log_masses = row_df["row_mass_log10"].values
    ones    = row_df["ones_count"].values

    # Map time → angle so the spiral winds through 1.5 full turns
    max_turns = 1.5
    angles_d = times / max(1, times[-1]) * max_turns * 2 * math.pi

    x_spiral = log_masses * __import__("numpy").cos(angles_d)
    y_spiral = log_masses * __import__("numpy").sin(angles_d)

    # Draw connecting spiral line in segments coloured by ones_count
    pts = __import__("numpy").array([x_spiral, y_spiral]).T.reshape(-1, 1, 2)
    segs = __import__("numpy").concatenate([pts[:-1], pts[1:]], axis=1)
    ones_norm = Normalize(vmin=int(ones.min()), vmax=int(ones.max()))
    lc = LineCollection(segs, cmap="cool", norm=ones_norm, lw=2.5, alpha=0.9)
    lc.set_array(ones[:-1])
    ax_spiral.add_collection(lc)
    ax_spiral.scatter(x_spiral, y_spiral, c=ones, cmap="cool", norm=ones_norm,
                      s=30, zorder=3, edgecolors="none")

    # Annotate first and last point
    ax_spiral.text(x_spiral[0], y_spiral[0], f" t=0", color="#88ffdd",
                   fontsize=8, va="bottom")
    ax_spiral.text(x_spiral[-1], y_spiral[-1], f" t={int(times[-1])} (halt)",
                   color="#ff88aa", fontsize=8, va="bottom")
    ax_spiral.axhline(0, color="#333344", lw=0.5)
    ax_spiral.axvline(0, color="#333344", lw=0.5)

    cb_d = fig.colorbar(lc, ax=ax_spiral, fraction=0.04, pad=0.04)
    cb_d.set_label("ones count on tape", color="#b0b0c0")
    plt.setp(cb_d.ax.yaxis.get_ticklabels(), color="#b0b0c0")
    cb_d.ax.yaxis.set_tick_params(color="#b0b0c0")

    ax_spiral.set_title(
        "D — Waring-Goldbach Row-Mass Spiral (log-polar)\n"
        "Angle ∝ time,  radius = log₁₀(N_t)  where  N_t = Σ_x p^k.\n"
        "The mass grows as computation accumulates → nautilus-shell spiral.",
        pad=6,
    )
    ax_spiral.set_xlabel("log₁₀(row mass) · cos(angle)")
    ax_spiral.set_ylabel("log₁₀(row mass) · sin(angle)")
    ax_spiral.set_aspect("equal")
    ax_spiral.margins(0.15)
    ax_spiral.grid(True, alpha=0.1, color="#444455")

    fig.suptitle(
        f"Busy-Beaver {machine.name.upper()} × Waring-Goldbach "
        "Space–Time Diagrams  |  Prime-Boost Shells\n"
        "Each (t, x) cell → prime p ≡ r (mod 30),  boost = p^k  where k encodes "
        "computational depth.  N_t = Σ_x p^k  (Waring-Goldbach row mass).",
        fontsize=12, color="#e8e8f0", y=1.01,
    )

    fig.savefig(out_path, dpi=180, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved shell diagram: {out_path}")


def build_output_paths(machine_name: str, args: argparse.Namespace) -> dict[str, str]:
    stem = f"busy_beaver_waring_goldbach_{machine_name}"
    return {
        "cells": args.cells_out or f"outputs/csv/{stem}_cells.csv",
        "rows": args.rows_out or f"outputs/csv/{stem}_rows.csv",
        "tiles": args.tiles_out or f"outputs/csv/{stem}_tiles.csv",
        "plot": args.plot_out or f"outputs/plots/{stem}_dashboard.png",
    }


def parse_args() -> argparse.Namespace:
    machines = sorted(builtin_machines())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--machine", choices=machines, default="bb3")
    parser.add_argument("--max-steps", type=int, default=64)
    parser.add_argument(
        "--padding",
        type=int,
        default=2,
        help="Extra tape cells shown on both sides of the active support.",
    )
    parser.add_argument("--cells-out", default=None)
    parser.add_argument("--rows-out", default=None)
    parser.add_argument("--tiles-out", default=None)
    parser.add_argument("--plot-out", default=None)
    parser.add_argument("--plot", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    machines = builtin_machines()
    machine = machines[args.machine]

    history = simulate_machine(machine, max_steps=args.max_steps)
    span = history_span(history, padding=max(1, args.padding))
    cell_records, cell_lookup = build_cell_records(history, machine, span)
    row_records = build_row_records(history, span, cell_lookup)
    tile_records = build_transition_plaquettes(history, machine, cell_lookup)

    cell_df = pd.DataFrame(cell_records)
    row_df = pd.DataFrame(row_records)
    tile_df = pd.DataFrame(tile_records)

    paths = build_output_paths(machine.name, args)
    for path in (paths["cells"], paths["rows"], paths["tiles"]):
        ensure_parent_dir(path)

    stringify_large_int_columns(cell_df, ["term_value"]).to_csv(paths["cells"], index=False)
    stringify_large_int_columns(row_df, ["row_mass"]).to_csv(paths["rows"], index=False)
    stringify_large_int_columns(tile_df, ["plaquette_power_sum"]).to_csv(paths["tiles"], index=False)

    if args.plot:
        plot_tableau(cell_df, row_df, machine, paths["plot"])

    halted = bool(history[-1].halted)
    step_count = len(history) - 1
    ones_count = sum(history[-1].tape.values())
    print(f"machine={machine.name}")
    print(f"halted={halted}")
    print(f"steps={step_count}")
    print(f"final_ones={ones_count}")
    print(f"visible_span=[{span[0]}, {span[-1]}]")
    print(f"cells_csv={paths['cells']}")
    print(f"rows_csv={paths['rows']}")
    print(f"tiles_csv={paths['tiles']}")
    if args.plot:
        print(f"plot_png={paths['plot']}")
    if not halted:
        print("warning=max_steps reached before halting")


if __name__ == "__main__":
    main()
