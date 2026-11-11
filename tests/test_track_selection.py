"""Independent examples distinguishing event multiplicity from one chosen track."""
import sys
import unittest
from dataclasses import replace
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from track_selection import Track, select
from analyze_thickness_mechanisms import classify_channel, select_agreed_track


class SelectionTest(unittest.TestCase):
    def setUp(self):
        self.a = Track(1, 15, 5., 0., 6., 140., .1, True, True)
        self.b = Track(2, 15, 2., 0., 2., 140., .1, True, True)

    def test_single(self):
        self.assertEqual(select([self.a], 'single'), self.a)
        self.assertIsNone(select([], 'single'))

    def test_extra_track_veto_before_cuts(self):
        self.assertIsNone(select([self.a, replace(self.b, sclose_cm=20.)], 'single'))

    def test_smallest_angle_is_not_largest_momentum(self):
        self.assertEqual(select([self.a, self.b], 'min_angle'), self.b)
        self.assertEqual(select([self.a, self.b], 'legacy'), self.a)
        self.assertEqual(select([replace(self.a, momentum_gev=.1),
                                 replace(self.b, momentum_gev=20.)], 'min_angle').track_id, 2)

    def test_do_not_search_for_passing_candidate(self):
        self.assertIsNone(select([self.a, replace(self.b, sclose_cm=20.)], 'min_angle'))

    def test_required_chambers_and_switch(self):
        self.assertIsNone(select([replace(self.a, chamber_mask=3)], 'min_angle'))
        self.assertIsNotNone(select([replace(self.a, chamber_mask=3, theta_deg=12.)], 'min_angle'))
        self.assertIsNone(select([self.a, replace(self.b, chamber_mask=3)], 'min_angle'))

    def test_strict_one_degree_eligibility_before_ranking(self):
        below = replace(self.a, theta_deg=.5)
        boundary = replace(self.a, track_id=3, theta_deg=1.)
        self.assertEqual(select([below, boundary, self.b], 'min_angle'), self.b)
        self.assertIsNone(select([below, boundary], 'min_angle'))
        self.assertEqual(select([below, self.b], 'min_angle_uncut'), below)

    def test_one_degree_does_not_change_multiplicity(self):
        self.assertIsNone(select([replace(self.a, theta_deg=.5), self.b], 'single'))
        self.assertEqual(select([replace(self.a, theta_deg=.5), self.b], 'min_angle'), self.b)

    def test_configurable_strict_threshold(self):
        self.assertEqual(select([self.a, self.b], 'min_angle', min_angle_deg=2.), self.a)

    def test_deterministic_angle_tie(self):
        self.assertEqual(select([replace(self.b, theta_deg=5.), self.a], 'min_angle'), self.a)

    def test_mechanism_channel_classification(self):
        primary = replace(self.a, root_track_id=1, generation=0)
        elastic_secondary = replace(
            self.a, track_id=2, root_track_id=1, generation=1,
            lineage_process_category=1)
        inelastic_secondary = replace(
            self.a, track_id=3, root_track_id=1, generation=2,
            lineage_process_category=2)
        self.assertEqual(classify_channel(primary, 1), 'primary_elastic')
        self.assertEqual(classify_channel(primary, 2), 'primary_inelastic')
        self.assertEqual(classify_channel(elastic_secondary, 0),
                         'secondary_elastic')
        self.assertEqual(classify_channel(inelastic_secondary, 0),
                         'secondary_inelastic')

    def test_thickness_policy_is_smallest_angle_without_momentum(self):
        high_momentum = replace(self.a, momentum_gev=100.)
        low_momentum = replace(self.b, momentum_gev=.01)
        self.assertEqual(
            select_agreed_track([high_momentum, low_momentum], 55.88),
            low_momentum,
        )

    def test_thickness_policy_filters_one_degree_before_ranking(self):
        below = replace(self.a, theta_deg=.5)
        boundary = replace(self.a, track_id=3, theta_deg=1.)
        self.assertEqual(select_agreed_track([below, boundary, self.b], 55.88),
                         self.b)

    def test_thickness_policy_always_requires_all_four_chambers(self):
        front_pair_only = replace(self.a, chamber_mask=3, theta_deg=20.)
        self.assertIsNone(select_agreed_track([front_pair_only], 55.88))

    def test_thickness_policy_does_not_replace_failed_smallest_track(self):
        failed_smallest = replace(self.b, chamber_mask=3)
        self.assertIsNone(select_agreed_track([self.a, failed_smallest], 55.88))


if __name__ == '__main__':
    unittest.main()
