# -*- coding: utf-8 -*-
"""Generate a 3D propeller + MRF CFD domain in Gmsh for OpenFOAM.

The geometry uses a NACA 9112 blade section and spanwise chord/twist control
points taken from the current propeller-design study in this repository.

Outputs:
- propeller_solid.step
- propeller_mrf_domain.msh

Physical groups:
- inlet
- outlet
- outerWall
- propeller
- rotatingZone
- stationaryZone
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

try:
    import gmsh
except ImportError as exc:  # pragma: no cover - import guard only
    raise SystemExit(
        "gmsh is required. Install it with `pip install gmsh` or use a Gmsh "
        "Python environment."
    ) from exc


# ============================================================
# 1. Propeller inputs
# ============================================================
R = 1.587
R_ROOT = 0.150
N_BLADES = 5

# Propeller axis is the x-axis. One blade is built along +z and then copied
# around the x-axis.
HUB_RADIUS = 0.155
HUB_LENGTH = 0.35
ROOT_BLEND_RADIUS = 0.11
ROOT_BLEND_CHORD_SCALE = 0.60

CONTROL_STATIONS = np.linspace(0.0, 1.0, 7)

DESIGN_PRESETS = {
    "ct": {
        "label": "max-CT",
        "chord_ctrl": np.linspace(0.40, 0.20, 7),
        "twist_ctrl_deg": np.array([87.3, 78.0, 69.5, 60.8, 53.8, 48.3, 44.0]),
    },
    "eta": {
        "label": "max-eta",
        "chord_ctrl": np.array([0.05, 0.045, 0.040, 0.035, 0.030, 0.025, 0.020]),
        "twist_ctrl_deg": np.array([86.0, 74.5, 64.0, 56.0, 49.8, 45.0, 38.5]),
    },
}

# Airfoil and meshing controls. These defaults bias toward smoother cell-size
# transitions, a less aggressive MRF interface, and a refined downstream wake
# while still staying in the range of an 8 GB VM-scale OpenFOAM test case.
AIRFOIL_CODE = "9112"
N_AIRFOIL_POINTS = 51
N_SPAN_SECTIONS = 8
PITCH_AXIS_FRACTION = 0.25

NEAR_BLADE_MESH_SIZE = 0.05
MID_DOMAIN_MESH_SIZE = 0.25
FARFIELD_MESH_SIZE = 0.75
TIP_MESH_SIZE = 0.05

NEAR_REFINEMENT_RADIUS = 0.10 * R
MID_REFINEMENT_RADIUS = 0.60 * R
INTERFACE_MESH_SIZE = 0.18
INTERFACE_REFINEMENT_THICKNESS = 0.12 * R
WAKE_MESH_SIZE = 0.18
WAKE_UPSTREAM_BUFFER = 0.10 * R
WAKE_DOWNSTREAM_LENGTH = 2.00 * R
WAKE_RADIUS = 0.70 * R

BL_N_LAYERS = 5
BL_FIRST_LAYER_THICKNESS = 0.001
BL_GROWTH_RATE = 1.2
BL_MIN_SURFACE_AREA = 5.0e-3
DEFAULT_USE_BOUNDARY_LAYER = False
DEFAULT_USE_EXPERIMENTAL_REFINEMENT = False

# CFD domain extents, all aligned with the x-axis for OpenFOAM MRF.
UPSTREAM_LENGTH = 3.0 * R
DOWNSTREAM_LENGTH = 6.0 * R
FARFIELD_RADIUS = 5.0 * R
MRF_RADIUS = 1.20 * R
MRF_HALF_LENGTH = 0.75 * R


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--design",
        choices=sorted(DESIGN_PRESETS),
        default="ct",
        help="Blade geometry preset to use.",
    )
    parser.add_argument(
        "--output-dir",
        default="gmsh_output",
        help="Directory for STEP and MSH outputs.",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Open the Gmsh GUI after building the mesh.",
    )
    parser.add_argument(
        "--experimental-refinement",
        action="store_true",
        default=DEFAULT_USE_EXPERIMENTAL_REFINEMENT,
        help="Enable additional wake and MRF-interface refinement fields.",
    )
    bl_group = parser.add_mutually_exclusive_group()
    bl_group.add_argument(
        "--boundary-layer",
        dest="use_boundary_layer",
        action="store_true",
        help="Enable prism-layer generation on the propeller surface.",
    )
    bl_group.add_argument(
        "--no-boundary-layer",
        dest="use_boundary_layer",
        action="store_false",
        help="Skip prism-layer generation and build a tetra-only volume mesh.",
    )
    parser.add_argument(
        "--strict-boundary-layer",
        action="store_true",
        help="Fail instead of falling back if prism-layer generation does not mesh cleanly.",
    )
    parser.set_defaults(use_boundary_layer=DEFAULT_USE_BOUNDARY_LAYER)
    return parser.parse_args()


def cosine_spacing(num_points: int) -> np.ndarray:
    beta = np.linspace(0.0, math.pi, num_points)
    return 0.5 * (1.0 - np.cos(beta))


def naca4_airfoil_surfaces(code: str, num_points: int) -> tuple[np.ndarray, np.ndarray]:
    if len(code) != 4 or not code.isdigit():
        raise ValueError(f"Expected a 4-digit NACA code, got {code!r}")

    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:]) / 100.0

    x = cosine_spacing(num_points)
    yt = 5.0 * t * (
        0.2969 * np.sqrt(np.maximum(x, 1e-12))
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1036 * x**4
    )

    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)

    before_p = x < p
    after_p = ~before_p

    yc[before_p] = m / p**2 * (2.0 * p * x[before_p] - x[before_p] ** 2)
    yc[after_p] = m / (1.0 - p) ** 2 * ((1.0 - 2.0 * p) + 2.0 * p * x[after_p] - x[after_p] ** 2)

    dyc_dx[before_p] = 2.0 * m / p**2 * (p - x[before_p])
    dyc_dx[after_p] = 2.0 * m / (1.0 - p) ** 2 * (p - x[after_p])

    theta = np.arctan(dyc_dx)

    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)

    upper = np.column_stack([xu[::-1], yu[::-1]])
    lower = np.column_stack([xl, yl])

    # Force exact closure at the leading and trailing edges so the loft can use
    # shared end points instead of a tiny trailing-edge line, which previously
    # produced overlapping-facet failures during 3D meshing.
    trailing_edge = np.array([1.0, 0.0])
    leading_edge = np.array([0.0, 0.0])
    upper[0] = trailing_edge
    lower[-1] = trailing_edge
    upper[-1] = leading_edge
    lower[0] = leading_edge
    return upper, lower


def interpolate_geometry(control_values: np.ndarray, span_fracs: np.ndarray) -> np.ndarray:
    return np.interp(span_fracs, CONTROL_STATIONS, control_values)


def build_section_points(
    radius: float,
    chord: float,
    twist_deg: float,
    airfoil_xy: np.ndarray,
) -> np.ndarray:
    x_local = (airfoil_xy[:, 0] - PITCH_AXIS_FRACTION) * chord
    y_local = airfoil_xy[:, 1] * chord

    beta = math.radians(twist_deg)
    x_rot = x_local * math.cos(beta) - y_local * math.sin(beta)
    y_rot = x_local * math.sin(beta) + y_local * math.cos(beta)
    z_rot = np.full_like(x_rot, radius)

    return np.column_stack([x_rot, y_rot, z_rot])


def add_airfoil_wire(
    upper_section_points: np.ndarray,
    lower_section_points: np.ndarray,
    mesh_size: float,
) -> int:
    trailing_edge = gmsh.model.occ.addPoint(
        float(upper_section_points[0, 0]),
        float(upper_section_points[0, 1]),
        float(upper_section_points[0, 2]),
        mesh_size,
    )
    upper_tags = [trailing_edge]
    upper_tags.extend(
        gmsh.model.occ.addPoint(float(x), float(y), float(z), mesh_size)
        for x, y, z in upper_section_points[1:-1]
    )

    leading_edge = gmsh.model.occ.addPoint(
        float(upper_section_points[-1, 0]),
        float(upper_section_points[-1, 1]),
        float(upper_section_points[-1, 2]),
        mesh_size,
    )
    upper_tags.append(leading_edge)

    lower_tags = [leading_edge]
    lower_tags.extend(
        gmsh.model.occ.addPoint(float(x), float(y), float(z), mesh_size)
        for x, y, z in lower_section_points[1:-1]
    )
    lower_tags.append(trailing_edge)

    upper_spline = gmsh.model.occ.addSpline(upper_tags)
    lower_spline = gmsh.model.occ.addSpline(lower_tags)
    return gmsh.model.occ.addWire([upper_spline, lower_spline], checkClosed=True)


def extract_dimtags(dimtags: list[tuple[int, int]], dim: int) -> list[tuple[int, int]]:
    return [entity for entity in dimtags if entity[0] == dim]


def boundary_surfaces(dimtags: list[tuple[int, int]]) -> list[int]:
    return sorted({tag for dim, tag in gmsh.model.getBoundary(dimtags, oriented=False, recursive=False) if dim == 2})


def shared_surfaces(
    dimtags_a: list[tuple[int, int]],
    dimtags_b: list[tuple[int, int]],
) -> list[int]:
    surfaces_a = set(boundary_surfaces(dimtags_a))
    surfaces_b = set(boundary_surfaces(dimtags_b))
    return sorted(surfaces_a & surfaces_b)


def entity_bbox(dim: int, tag: int) -> tuple[float, float, float, float, float, float]:
    return gmsh.model.getBoundingBox(dim, tag)


def surface_bbox(tag: int) -> tuple[float, float, float, float, float, float]:
    return entity_bbox(2, tag)


def is_plane_at_x(tag: int, x_target: float, tol: float = 1e-5) -> bool:
    xmin, _, _, xmax, _, _ = surface_bbox(tag)
    return abs(xmin - x_target) < tol and abs(xmax - x_target) < tol


def max_transverse_extent_from_surface_bbox(tag: int) -> float:
    _, ymin, zmin, _, ymax, zmax = surface_bbox(tag)
    return max(abs(ymin), abs(ymax), abs(zmin), abs(zmax))


def max_transverse_extent_from_entity_bbox(dim: int, tag: int) -> float:
    _, ymin, zmin, _, ymax, zmax = entity_bbox(dim, tag)
    return max(abs(ymin), abs(ymax), abs(zmin), abs(zmax))


def touches_x_plane(x_value: float, target: float, tol: float = 1e-5) -> bool:
    return abs(x_value - target) < tol


def cumulative_layer_heights(first_layer: float, growth_rate: float, num_layers: int) -> list[float]:
    heights = []
    total = 0.0
    thickness = first_layer
    for _ in range(num_layers):
        total += thickness
        heights.append(total)
        thickness *= growth_rate
    return heights


def configure_mesh_options(use_boundary_layer: bool) -> None:
    min_size = BL_FIRST_LAYER_THICKNESS if use_boundary_layer else min(NEAR_BLADE_MESH_SIZE, TIP_MESH_SIZE)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", FARFIELD_MESH_SIZE)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.Smoothing", 10)


def set_mesh_fields(
    propeller_surfaces: list[int],
    interface_surfaces: list[int],
    use_experimental_refinement: bool,
) -> None:
    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "FacesList", propeller_surfaces)

    near_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(near_field, "InField", distance)
    gmsh.model.mesh.field.setNumber(near_field, "SizeMin", NEAR_BLADE_MESH_SIZE)
    gmsh.model.mesh.field.setNumber(near_field, "SizeMax", MID_DOMAIN_MESH_SIZE)
    gmsh.model.mesh.field.setNumber(near_field, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(near_field, "DistMax", NEAR_REFINEMENT_RADIUS)

    mid_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(mid_field, "InField", distance)
    gmsh.model.mesh.field.setNumber(mid_field, "SizeMin", MID_DOMAIN_MESH_SIZE)
    gmsh.model.mesh.field.setNumber(mid_field, "SizeMax", FARFIELD_MESH_SIZE)
    gmsh.model.mesh.field.setNumber(mid_field, "DistMin", NEAR_REFINEMENT_RADIUS)
    gmsh.model.mesh.field.setNumber(mid_field, "DistMax", MID_REFINEMENT_RADIUS)

    fields = [near_field, mid_field]

    if use_experimental_refinement and interface_surfaces:
        interface_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(interface_distance, "FacesList", interface_surfaces)

        interface_field = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(interface_field, "InField", interface_distance)
        gmsh.model.mesh.field.setNumber(interface_field, "SizeMin", INTERFACE_MESH_SIZE)
        gmsh.model.mesh.field.setNumber(interface_field, "SizeMax", MID_DOMAIN_MESH_SIZE)
        gmsh.model.mesh.field.setNumber(interface_field, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(interface_field, "DistMax", INTERFACE_REFINEMENT_THICKNESS)
        fields.append(interface_field)

    if use_experimental_refinement:
        wake_field = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(wake_field, "VIn", WAKE_MESH_SIZE)
        gmsh.model.mesh.field.setNumber(wake_field, "VOut", FARFIELD_MESH_SIZE)
        gmsh.model.mesh.field.setNumber(wake_field, "XMin", -WAKE_UPSTREAM_BUFFER)
        gmsh.model.mesh.field.setNumber(wake_field, "XMax", WAKE_DOWNSTREAM_LENGTH)
        gmsh.model.mesh.field.setNumber(wake_field, "YMin", -WAKE_RADIUS)
        gmsh.model.mesh.field.setNumber(wake_field, "YMax", WAKE_RADIUS)
        gmsh.model.mesh.field.setNumber(wake_field, "ZMin", -WAKE_RADIUS)
        gmsh.model.mesh.field.setNumber(wake_field, "ZMax", WAKE_RADIUS)
        fields.append(wake_field)

    background = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(background, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(background)


def add_named_physical_group(dim: int, tags: list[int], name: str) -> None:
    if not tags:
        raise RuntimeError(f"Physical group {name!r} would be empty.")
    phys = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, phys, name)


def optimize_volume_mesh() -> None:
    gmsh.model.mesh.optimize()
    gmsh.model.mesh.optimize("Netgen")


def build_propeller_geometry(design: dict[str, np.ndarray]) -> list[tuple[int, int]]:
    upper_airfoil_xy, lower_airfoil_xy = naca4_airfoil_surfaces(AIRFOIL_CODE, N_AIRFOIL_POINTS)

    section_fracs = np.linspace(0.0, 1.0, N_SPAN_SECTIONS)
    section_radii = R_ROOT + section_fracs * (R - R_ROOT)
    section_chords = interpolate_geometry(np.asarray(design["chord_ctrl"]), section_fracs)
    section_twists = interpolate_geometry(np.asarray(design["twist_ctrl_deg"]), section_fracs)

    # Add one extra inboard section inside the hub to create a stronger blade/hub
    # overlap. Without this blend section, the root intersection can be too sharp
    # for the 3D mesher and cause PLC recovery failures.
    section_radii = np.concatenate([[ROOT_BLEND_RADIUS], section_radii])
    section_chords = np.concatenate([[ROOT_BLEND_CHORD_SCALE * section_chords[0]], section_chords])
    section_twists = np.concatenate([[section_twists[0]], section_twists])

    wire_tags = []
    root_to_tip_fraction = (section_radii - section_radii.min()) / (section_radii.max() - section_radii.min())
    for frac, radius, chord, twist in zip(root_to_tip_fraction, section_radii, section_chords, section_twists):
        mesh_size = NEAR_BLADE_MESH_SIZE if frac < 0.9 else TIP_MESH_SIZE
        upper_section_points = build_section_points(radius, chord, twist, upper_airfoil_xy)
        lower_section_points = build_section_points(radius, chord, twist, lower_airfoil_xy)
        wire_tags.append(add_airfoil_wire(upper_section_points, lower_section_points, mesh_size))

    loft_entities = gmsh.model.occ.addThruSections(wire_tags, makeSolid=True, makeRuled=False)
    base_blade = extract_dimtags(loft_entities, 3)
    if not base_blade:
        raise RuntimeError("Failed to build blade loft volume from section wires.")

    all_blades = list(base_blade)
    for blade_idx in range(1, N_BLADES):
        copied = gmsh.model.occ.copy(base_blade)
        gmsh.model.occ.rotate(copied, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, blade_idx * 2.0 * math.pi / N_BLADES)
        all_blades.extend(copied)

    hub = [(3, gmsh.model.occ.addCylinder(-0.5 * HUB_LENGTH, 0.0, 0.0, HUB_LENGTH, 0.0, 0.0, HUB_RADIUS))]
    propeller_entities, _ = gmsh.model.occ.fuse(hub, all_blades, removeObject=True, removeTool=True)
    propeller_volumes = extract_dimtags(propeller_entities, 3)
    if not propeller_volumes:
        raise RuntimeError("Failed to fuse hub and blades into a propeller solid.")

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    # OCC booleans and duplicate cleanup can renumber or replace the fused
    # solid entities. Refresh the live 3D tags before returning them to the
    # downstream cut/fragment operations.
    propeller_volumes = gmsh.model.getEntities(3)
    if not propeller_volumes:
        raise RuntimeError("Propeller solid disappeared after OCC duplicate cleanup.")

    return propeller_volumes


def build_mrf_domain(
    propeller_volumes: list[tuple[int, int]],
) -> tuple[list[tuple[int, int]], list[tuple[int, int]], list[int], list[int], list[int], list[int], list[int]]:
    farfield_tag = gmsh.model.occ.addCylinder(-UPSTREAM_LENGTH, 0.0, 0.0, UPSTREAM_LENGTH + DOWNSTREAM_LENGTH, 0.0, 0.0, FARFIELD_RADIUS)
    mrf_tag = gmsh.model.occ.addCylinder(-MRF_HALF_LENGTH, 0.0, 0.0, 2.0 * MRF_HALF_LENGTH, 0.0, 0.0, MRF_RADIUS)

    # Remove the propeller solid from the rotating fluid zone. The propeller
    # itself should not remain as a meshed volume for OpenFOAM MRF; we only
    # want the blade/hub cavity surfaces to become wall patches.
    rotating_zone, _ = gmsh.model.occ.cut([(3, mrf_tag)], propeller_volumes, removeObject=True, removeTool=True)
    stationary_zone, _ = gmsh.model.occ.cut([(3, farfield_tag)], rotating_zone, removeObject=True, removeTool=False)

    rotating_volumes = extract_dimtags(rotating_zone, 3)
    stationary_volumes = extract_dimtags(stationary_zone, 3)
    if not rotating_volumes or not stationary_volumes:
        raise RuntimeError("Failed to build rotating/stationary CFD zones.")

    # Make the rotating and stationary fluid zones conformal at their mutual
    # interface. Without this fragment, Gmsh can keep coincident but distinct
    # interface surfaces, which later show up as overlapping facets in 3D
    # meshing.
    gmsh.model.occ.fragment(rotating_volumes, stationary_volumes, removeObject=True, removeTool=True)
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    residual_solids = []
    rotating_volumes = []
    stationary_volumes = []
    for dim, tag in gmsh.model.getEntities(3):
        xmin, _, _, xmax, _, _ = entity_bbox(dim, tag)
        max_radius = max_transverse_extent_from_entity_bbox(dim, tag)
        inside_mrf_bbox = (
            xmin >= -MRF_HALF_LENGTH - 1e-5
            and xmax <= MRF_HALF_LENGTH + 1e-5
            and max_radius <= MRF_RADIUS + 1e-5
        )
        if inside_mrf_bbox:
            # The actual rotating fluid zone reaches both axial MRF end planes and
            # the cylindrical shell. Smaller volumes that remain inside the MRF
            # box are stray OCC solids from the blade/hub booleans and must not be
            # passed to the mesher.
            spans_mrf_length = (
                touches_x_plane(xmin, -MRF_HALF_LENGTH) and touches_x_plane(xmax, MRF_HALF_LENGTH)
            )
            reaches_mrf_shell = max_radius >= 0.95 * MRF_RADIUS
            if spans_mrf_length or reaches_mrf_shell:
                rotating_volumes.append((dim, tag))
            else:
                residual_solids.append((dim, tag))
        else:
            stationary_volumes.append((dim, tag))

    if residual_solids:
        gmsh.model.occ.remove(residual_solids, recursive=True)
        gmsh.model.occ.synchronize()
        rotating_volumes = [entity for entity in rotating_volumes if entity not in residual_solids]

    if not rotating_volumes or not stationary_volumes:
        raise RuntimeError("Failed to classify rotating/stationary CFD volumes after boolean operations.")

    rotating_surfaces = boundary_surfaces(rotating_volumes)
    stationary_surfaces = boundary_surfaces(stationary_volumes)
    interface_surfaces = sorted(shared_surfaces(rotating_volumes, stationary_volumes))
    interface_surface_set = set(interface_surfaces)
    inlet, outlet, outer_wall = [], [], []
    for surface_tag in stationary_surfaces:
        if is_plane_at_x(surface_tag, -UPSTREAM_LENGTH):
            inlet.append(surface_tag)
        elif is_plane_at_x(surface_tag, DOWNSTREAM_LENGTH):
            outlet.append(surface_tag)
        elif max_transverse_extent_from_surface_bbox(surface_tag) > 0.95 * FARFIELD_RADIUS:
            outer_wall.append(surface_tag)

    propeller_surfaces = []
    for surface_tag in rotating_surfaces:
        if is_plane_at_x(surface_tag, -MRF_HALF_LENGTH):
            continue
        if is_plane_at_x(surface_tag, MRF_HALF_LENGTH):
            continue
        if surface_tag in interface_surface_set:
            continue
        if max_transverse_extent_from_surface_bbox(surface_tag) > 0.95 * MRF_RADIUS:
            continue
        propeller_surfaces.append(surface_tag)

    if not propeller_surfaces:
        excluded_surfaces = interface_surface_set | set(inlet) | set(outlet) | set(outer_wall)
        for _, surface_tag in gmsh.model.getEntities(2):
            if surface_tag in excluded_surfaces:
                continue
            if is_plane_at_x(surface_tag, -MRF_HALF_LENGTH) or is_plane_at_x(surface_tag, MRF_HALF_LENGTH):
                continue
            if max_transverse_extent_from_surface_bbox(surface_tag) > 0.95 * FARFIELD_RADIUS:
                continue
            propeller_surfaces.append(surface_tag)

    return (
        rotating_volumes,
        stationary_volumes,
        interface_surfaces,
        propeller_surfaces,
        inlet,
        outlet,
        outer_wall,
    )


def should_extrude_boundary_layer(surface_tag: int) -> bool:
    if gmsh.model.getType(2, surface_tag) == "Plane":
        return False
    return gmsh.model.occ.getMass(2, surface_tag) >= BL_MIN_SURFACE_AREA


def add_propeller_boundary_layer(
    rotating_volumes: list[tuple[int, int]],
    propeller_surfaces: list[int],
) -> list[tuple[int, int]]:
    rotating_shell_surfaces = [tag for tag in boundary_surfaces(rotating_volumes) if tag not in propeller_surfaces]
    bl_source_surfaces = [tag for tag in propeller_surfaces if should_extrude_boundary_layer(tag)]
    skipped_surfaces = [tag for tag in propeller_surfaces if tag not in bl_source_surfaces]
    if not bl_source_surfaces:
        raise RuntimeError("No propeller surfaces qualified for prism-layer extrusion.")

    gmsh.model.occ.remove(rotating_volumes, recursive=False)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Geometry.ExtrudeReturnLateralEntities", 0)

    extruded = gmsh.model.geo.extrudeBoundaryLayer(
        [(2, tag) for tag in bl_source_surfaces],
        [1] * BL_N_LAYERS,
        cumulative_layer_heights(BL_FIRST_LAYER_THICKNESS, BL_GROWTH_RATE, BL_N_LAYERS),
        True,
    )

    top_surfaces = []
    boundary_layer_volumes = []
    for idx in range(1, len(extruded)):
        if extruded[idx][0] == 3:
            boundary_layer_volumes.append((3, extruded[idx][1]))
            if extruded[idx - 1][0] == 2:
                top_surfaces.append(extruded[idx - 1][1])

    gmsh.model.geo.synchronize()

    core_surface_loop = gmsh.model.geo.addSurfaceLoop(rotating_shell_surfaces + skipped_surfaces + top_surfaces)
    core_volume = gmsh.model.geo.addVolume([core_surface_loop])
    gmsh.model.geo.synchronize()

    return [(3, core_volume), *boundary_layer_volumes]


def build_and_mesh_case(
    design_key: str,
    output_dir: Path,
    use_boundary_layer: bool,
    use_experimental_refinement: bool,
    open_gui: bool,
) -> None:
    design = DESIGN_PRESETS[design_key]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("propeller_mrf")

    try:
        propeller_volumes = build_propeller_geometry(design)
        gmsh.write(str(output_dir / "propeller_solid.step"))

        (
            rotating_volumes,
            stationary_volumes,
            interface_surfaces,
            propeller_surfaces,
            inlet,
            outlet,
            outer_wall,
        ) = build_mrf_domain(propeller_volumes)
        if use_boundary_layer:
            rotating_volumes = add_propeller_boundary_layer(rotating_volumes, propeller_surfaces)

        add_named_physical_group(3, [tag for _, tag in rotating_volumes], "rotatingZone")
        add_named_physical_group(3, [tag for _, tag in stationary_volumes], "stationaryZone")

        add_named_physical_group(2, propeller_surfaces, "propeller")
        add_named_physical_group(2, inlet, "inlet")
        add_named_physical_group(2, outlet, "outlet")
        add_named_physical_group(2, outer_wall, "outerWall")

        configure_mesh_options(use_boundary_layer)
        set_mesh_fields(propeller_surfaces, interface_surfaces, use_experimental_refinement)
        gmsh.model.mesh.generate(3)
        optimize_volume_mesh()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(str(output_dir / "propeller_mrf_domain.msh"))

        if open_gui:
            gmsh.fltk.run()
    finally:
        gmsh.finalize()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    use_boundary_layer = args.use_boundary_layer

    try:
        build_and_mesh_case(
            args.design,
            output_dir,
            use_boundary_layer,
            args.experimental_refinement,
            args.gui,
        )
    except Exception as exc:
        if not use_boundary_layer or args.strict_boundary_layer:
            raise
        print(
            "Warning: prism-layer meshing failed in Gmsh; retrying with the same "
            "refinement zones but without boundary-layer extrusion.",
            flush=True,
        )
        build_and_mesh_case(
            args.design,
            output_dir,
            False,
            args.experimental_refinement,
            args.gui,
        )


if __name__ == "__main__":
    main()
