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
