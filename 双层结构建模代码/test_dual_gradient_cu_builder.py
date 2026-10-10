import importlib.util
import json
import tempfile
import sys
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("dual_gradient_cu_builder.py")
SPEC = importlib.util.spec_from_file_location("dual_gradient_cu_builder", MODULE_PATH)
builder = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = builder
SPEC.loader.exec_module(builder)


class NanotwinPlanningTests(unittest.TestCase):
    def hcp_positions(self, labels):
        n = len(labels)
        return [i for i in range(n) if labels[(i - 1) % n] == labels[(i + 1) % n]]

    def assert_uniform_periodic_hcp_spacing(self, labels):
        hcp = self.hcp_positions(labels)
        self.assertGreaterEqual(len(hcp), 2)
        spacings = [hcp[i + 1] - hcp[i] for i in range(len(hcp) - 1)]
        wrap = (hcp[0] + len(labels)) - hcp[-1]
        self.assertEqual(len(set(spacings + [wrap])), 1)

    def test_stack_plan_can_start_with_twin_boundary_at_first_layer(self):
        plan = builder.plan_nanotwin_stack(
            total_layers=115,
            layers_per_twin=10,
        )

        self.assertEqual(plan["bottom_margin"], 0)
        self.assertEqual(plan["top_margin"], 5)
        self.assertEqual(plan["boundary_positions"][0], 0)

    def test_stack_plan_keeps_uniform_internal_spacing_and_puts_remainder_at_top(self):
        plan = builder.plan_nanotwin_stack(
            total_layers=115,
            layers_per_twin=10,
        )

        self.assertTrue(plan["boundary_positions"])
        self.assertEqual(
            sorted(set(plan["boundary_spacings"])),
            [10],
        )
        self.assertEqual(plan["boundary_positions"][-1], 110)

    def test_blocks_keep_their_own_start_phase_instead_of_continuing_previous_block(self):
        first = builder.BlockSpec(
            grain_code="G3x2",
            twin_code="T10",
            thickness_nm=12.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.5,
        )
        second = builder.BlockSpec(
            grain_code="G3x2",
            twin_code="T15",
            thickness_nm=10.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.5,
        )

        adjusted = builder.attach_block_edge_phases([first, second])

        self.assertEqual(adjusted[0].start_label, 0)
        self.assertEqual(adjusted[0].start_direction, 1)
        self.assertEqual(adjusted[1].start_label, 0)
        self.assertEqual(adjusted[1].start_direction, 1)

    def test_ttg_family_continues_twin_phase_between_blocks(self):
        block_a = builder.BlockSpec(
            grain_code="single",
            twin_code="T6",
            thickness_nm=24.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )
        block_b = builder.BlockSpec(
            grain_code="single",
            twin_code="T24",
            thickness_nm=20.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )

        adjusted = builder.attach_family_specific_phase_rules([block_a, block_b], family="TTG", a0=3.615)
        self.assertEqual(adjusted[0].start_label, 0)
        self.assertEqual(adjusted[0].start_direction, 1)
        self.assertEqual(
            (adjusted[1].start_label, adjusted[1].start_direction),
            (1, -1),
        )

    def test_ttg_family_uses_interface_spacing_search_for_t4_t30_ab(self):
        block_a = builder.BlockSpec(
            grain_code="single",
            twin_code="T4",
            thickness_nm=24.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )
        block_b = builder.BlockSpec(
            grain_code="single",
            twin_code="T30",
            thickness_nm=20.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )

        adjusted = builder.attach_family_specific_phase_rules([block_a, block_b], family="TTG", a0=3.615)

        self.assertEqual((adjusted[1].start_label, adjusted[1].start_direction), (2, -1))

    def test_ttg_family_uses_interface_spacing_search_for_t10_t15_ab(self):
        block_a = builder.BlockSpec(
            grain_code="single",
            twin_code="T10",
            thickness_nm=24.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )
        block_b = builder.BlockSpec(
            grain_code="single",
            twin_code="T15",
            thickness_nm=20.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )

        adjusted = builder.attach_family_specific_phase_rules([block_a, block_b], family="TTG", a0=3.615)

        self.assertEqual((adjusted[1].start_label, adjusted[1].start_direction), (0, 1))

    def test_selected_ttg_ab_phase_makes_interface_gap_match_a_or_b_period(self):
        block_a = builder.BlockSpec(
            grain_code="single",
            twin_code="T10",
            thickness_nm=24.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )
        block_b = builder.BlockSpec(
            grain_code="single",
            twin_code="T15",
            thickness_nm=20.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )

        selected_label, selected_direction = builder.select_ttg_ab_interface_phase(block_a, block_b, 3.615)
        d111 = 3.615 / (3 ** 0.5)
        total_a = builder.choose_compatible_total_layers(
            max(2, int(round(builder.nm_to_angstrom(24.0) / d111))),
            builder.TWIN_N111["T10"],
        )
        total_b = builder.choose_compatible_total_layers(
            max(2, int(round(builder.nm_to_angstrom(20.0) / d111))),
            builder.TWIN_N111["T15"],
        )
        plan_a = builder.plan_nanotwin_stack(total_a, builder.TWIN_N111["T10"], 0, 1)
        plan_b = builder.plan_nanotwin_stack(total_b, builder.TWIN_N111["T15"], selected_label, selected_direction)
        gap = builder.interface_hcp_gap_layers(list(plan_a["labels"]), list(plan_b["labels"]))

        self.assertIn(gap, {builder.TWIN_N111["T10"], builder.TWIN_N111["T15"]})

    def test_ttg_abab_phase_search_keeps_all_interfaces_legal_for_ttg1(self):
        specs = [
            builder.BlockSpec("single", "T10", 12.0, 18.0, 1.0, 0.25, 0.35),
            builder.BlockSpec("single", "T15", 10.0, 18.0, 1.0, 0.25, 0.35),
            builder.BlockSpec("single", "T10", 12.0, 18.0, 1.0, 0.25, 0.35),
            builder.BlockSpec("single", "T15", 10.0, 18.0, 1.0, 0.25, 0.35),
        ]

        adjusted = builder.attach_family_specific_phase_rules(specs, family="TTG", a0=3.615)

        def labels_for(spec):
            d111 = 3.615 / (3 ** 0.5)
            total = builder.choose_compatible_total_layers(
                max(2, int(round(builder.nm_to_angstrom(spec.thickness_nm) / d111))),
                builder.TWIN_N111[spec.twin_code],
            )
            plan = builder.plan_nanotwin_stack(
                total,
                builder.TWIN_N111[spec.twin_code],
                spec.start_label,
                spec.start_direction,
            )
            return list(plan["labels"])

        labels = [labels_for(spec) for spec in adjusted]
        gap_1 = builder.interface_hcp_gap_layers(labels[0], labels[1])
        gap_2 = builder.interface_hcp_gap_layers(labels[1], labels[2])
        gap_3 = builder.interface_hcp_gap_layers(labels[2], labels[3])

        self.assertIn(gap_1, {10, 15})
        self.assertIn(gap_2, {10, 15})
        self.assertIn(gap_3, {10, 15})

    def test_short_block_uses_top_remainder_only(self):
        plan = builder.plan_nanotwin_stack(
            total_layers=57,
            layers_per_twin=10,
        )

        self.assertEqual(plan["bottom_margin"], 0)
        self.assertEqual(plan["top_margin"], 7)
        self.assertEqual(plan["boundary_positions"], [0, 10, 20, 30, 40, 50])

    def test_compatible_total_layers_snap_to_periodic_value(self):
        self.assertEqual(builder.choose_compatible_total_layers(115, 10), 120)
        self.assertEqual(builder.choose_compatible_total_layers(96, 15), 90)
        self.assertEqual(builder.choose_compatible_total_layers(57, 10), 60)
        self.assertEqual(builder.choose_compatible_total_layers(48, 15), 45)

    def test_snapped_total_layers_give_uniform_periodic_hcp_spacing(self):
        total_layers = builder.choose_compatible_total_layers(115, 10)
        plan = builder.plan_nanotwin_stack(total_layers=total_layers, layers_per_twin=10)

        self.assert_uniform_periodic_hcp_spacing(plan["labels"])

    def test_parser_uses_shifted_default_y_seed_offset(self):
        parser = builder.build_parser()
        args = parser.parse_args(
            ["block", "--grain", "G3x2", "--twin", "T10", "--thickness", "24"]
        )

        self.assertEqual(args.shifty, 0.35)

    def test_ttq_variant_rotates_b_layers_by_90_deg(self):
        model = {
            "grain_a": "G3x2",
            "grain_b": "G3x2",
            "twin_a": "T10",
            "twin_b": "T15",
            "family": "TTG",
        }
        context = builder.BuildContext(atomsk="dummy", root_dir=Path("."))

        specs = builder.gradient_variant_specs(model, "ABAB", context)

        self.assertEqual([spec.grain_code for spec in specs], ["single", "single", "single", "single"])
        self.assertEqual([spec.twin_code for spec in specs], ["T10", "T15", "T10", "T15"])
        self.assertEqual([spec.pattern_rotation_deg for spec in specs], [0.0, 0.0, 0.0, 0.0])
        self.assertEqual([spec.grain_angle_offset_deg for spec in specs], [0.0, 0.0, 0.0, 0.0])

    def test_gsg_variant_is_pure_grain_without_twins(self):
        model = {
            "grain_a": "G4x2",
            "grain_b": "G2x2",
            "twin_a": "T15",
            "twin_b": "T15",
            "family": "GSG",
        }
        context = builder.BuildContext(atomsk="dummy", root_dir=Path("."))

        specs = builder.gradient_variant_specs(model, "AB", context)

        self.assertEqual([spec.grain_code for spec in specs], ["G4x2", "G2x2"])
        self.assertEqual([spec.twin_code for spec in specs], [None, None])
        self.assertEqual([spec.pattern_rotation_deg for spec in specs], [0.0, 90.0])
        self.assertEqual([spec.grain_angle_offset_deg for spec in specs], [0.0, 90.0])

    def test_dg_variant_keeps_both_grains_and_twins(self):
        model = {
            "grain_a": "G4x2",
            "grain_b": "G2x2",
            "twin_a": "T10",
            "twin_b": "T20",
            "family": "DG",
        }
        context = builder.BuildContext(atomsk="dummy", root_dir=Path("."))

        specs = builder.gradient_variant_specs(model, "AB", context)

        self.assertEqual([spec.grain_code for spec in specs], ["G4x2", "G2x2"])
        self.assertEqual([spec.twin_code for spec in specs], ["T10", "T20"])
        self.assertEqual([spec.pattern_rotation_deg for spec in specs], [0.0, 90.0])
        self.assertEqual([spec.grain_angle_offset_deg for spec in specs], [0.0, 90.0])

    def test_collect_pure_components_includes_pure_twin_and_pure_grain(self):
        context = builder.BuildContext(atomsk="dummy", root_dir=Path("."))
        pairs = builder.collect_pure_components(builder.PLAN_MODELS, context)

        self.assertIn(("single", "T10"), pairs)
        self.assertIn(("G4x2", None), pairs)

    def test_classify_component_separates_pure_and_mixed_cases(self):
        self.assertEqual(
            builder.classify_component("single", "T15"),
            ("PURE_TWIN", "single_T15", "PURE_TWIN__single__T15"),
        )
        self.assertEqual(
            builder.classify_component("G4x2", None),
            ("PURE_GRAIN", "G4x2_FCC", "PURE_GRAIN__G4x2__FCC"),
        )
        self.assertEqual(
            builder.classify_component("G4x2", "T10"),
            ("COMPONENT", "G4x2_T10", "COMPONENT__G4x2__T10"),
        )

    def test_single_twin_block_bypasses_polycrystal_path(self):
        pure_twin = builder.BlockSpec(
            grain_code="single",
            twin_code="T15",
            thickness_nm=44.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )
        poly_twin = builder.BlockSpec(
            grain_code="G3x2",
            twin_code="T15",
            thickness_nm=20.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )

        self.assertTrue(builder.should_use_direct_block_build(pure_twin))
        self.assertFalse(builder.should_use_direct_block_build(poly_twin))

    def test_single_block_model_bypasses_merge_path(self):
        self.assertTrue(builder.should_copy_single_block_model(1))
        self.assertFalse(builder.should_copy_single_block_model(2))

    def test_single_twin_block_requires_unskew_for_lammps_export(self):
        pure_twin = builder.BlockSpec(
            grain_code="single",
            twin_code="T15",
            thickness_nm=44.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )
        pure_grain = builder.BlockSpec(
            grain_code="G4x2",
            twin_code=None,
            thickness_nm=44.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )

        self.assertTrue(builder.should_unskew_for_lammps_export(pure_twin))
        self.assertFalse(builder.should_unskew_for_lammps_export(pure_grain))

    def test_orthogonal_single_twin_builder_chooses_near_target_inplane_size(self):
        nx, ny = builder.choose_orthogonal_single_twin_repeats(a0=3.615, box_xy_nm=18.0)

        self.assertEqual((nx, ny), (41, 70))

    def test_parser_accepts_variant_filter_for_plan(self):
        parser = builder.build_parser()
        args = parser.parse_args(["plan", "--phase", "P0", "--models", "TTG-1", "--variants", "AB"])

        self.assertEqual(args.variants, "AB")

    def test_merge_paths_use_absolute_inputs_and_local_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            model_dir = root / "models" / "TTG" / "TTG-1" / "AB"
            block_dirs = [
                root / "cache" / "blocks" / "a",
                root / "cache" / "blocks" / "b",
            ]
            for block_dir in block_dirs:
                block_dir.mkdir(parents=True, exist_ok=True)

            inputs, merged_name, wrapped_name = builder.prepare_merge_paths(block_dirs, model_dir)

            self.assertTrue(all(not Path(item).is_absolute() for item in inputs))
            self.assertEqual(inputs[0], str(Path("..") / ".." / ".." / ".." / "cache" / "blocks" / "a" / "wrapped.cfg"))
            self.assertEqual(merged_name, "merged.cfg")
            self.assertEqual(wrapped_name, "wrapped.cfg")

    def test_single_twin_block_cache_requires_orthogonal_metadata(self):
        spec = builder.BlockSpec(
            grain_code="single",
            twin_code="T20",
            thickness_nm=20.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )

        with tempfile.TemporaryDirectory() as tmp:
            block_dir = Path(tmp)
            (block_dir / "wrapped.lmp").write_text("placeholder", encoding="utf-8")
            (block_dir / "metadata.json").write_text(
                json.dumps({"seed": {"mode": "nanotwin"}}),
                encoding="utf-8",
            )

            self.assertFalse(builder.block_cache_is_current(block_dir, spec))

            (block_dir / "metadata.json").write_text(
                json.dumps(
                    {
                        "seed": {
                            "mode": "orthogonal_nanotwin",
                            "orthogonal_single_twin": True,
                            "start_label": 1,
                            "start_direction": 1,
                        }
                    }
                ),
                encoding="utf-8",
            )
            self.assertFalse(builder.block_cache_is_current(block_dir, spec))

            (block_dir / "metadata.json").write_text(
                json.dumps(
                    {
                        "seed": {
                            "mode": "orthogonal_nanotwin",
                            "orthogonal_single_twin": True,
                            "start_label": 0,
                            "start_direction": 1,
                        }
                    }
                ),
                encoding="utf-8",
            )
            self.assertTrue(builder.block_cache_is_current(block_dir, spec))

    def test_orthogonal_single_twin_seed_preserves_requested_phase(self):
        with tempfile.TemporaryDirectory() as tmp:
            seed_path = Path(tmp) / "seed.xsf"
            info = builder.write_orthogonal_single_twin_seed(
                output_path=seed_path,
                a0=3.615,
                element="Cu",
                total_thickness_angstrom=200.0,
                twin_code="T30",
                box_xy_nm=18.0,
                start_label=1,
                start_direction=-1,
            )

        self.assertEqual(info["start_label"], 1)
        self.assertEqual(info["start_direction"], -1)

    def test_model_output_is_invalid_when_cached_block_ids_do_not_match(self):
        spec = builder.BlockSpec(
            grain_code="single",
            twin_code="T10",
            thickness_nm=24.0,
            box_xy_nm=18.0,
            packing_scale=1.0,
            shiftx=0.25,
            shifty=0.35,
        )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            model_dir = root / "models" / "TTG" / "TTG-1" / "AB"
            model_dir.mkdir(parents=True, exist_ok=True)
            (model_dir / "wrapped.lmp").write_text("placeholder", encoding="utf-8")
            (model_dir / "metadata.json").write_text(
                json.dumps(
                    {
                        "family": "TTG",
                        "phase": "P0",
                        "variant": "AB",
                        "layers": [builder.asdict(spec)],
                        "block_dirs": ["cache\\blocks\\stale-block-id"],
                        "single_block_copy": True,
                    }
                ),
                encoding="utf-8",
            )

            self.assertFalse(
                builder.model_output_is_current(
                    model_dir=model_dir,
                    effective_specs=[spec],
                    family="TTG",
                    phase="P0",
                    variant="AB",
                    root_dir=root,
                )
            )


if __name__ == "__main__":
    unittest.main()
