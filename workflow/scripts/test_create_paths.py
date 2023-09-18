import unittest, importlib
import tempfile
import io
import gzip
from contextlib import redirect_stdout
from create_paths import Variant, Cluster, variants_overlap, group_variants, print_vcf, run_paths


class Test_Variant(unittest.TestCase):
	def test_variant(self):
		variant = Variant("chr1", 100, "A", "G", 0)
		self.assertEqual(variant.chrom(), "chr1")
		self.assertEqual(variant.start(), 100)
		self.assertEqual(variant.ref(), "A")
		self.assertEqual(variant.end(), 101)
		self.assertEqual(variant.alt(), "G")
		self.assertEqual(variant.index(), 0)
		self.assertEqual(variant.id(), '.')


class Test_variants_overlap(unittest.TestCase):
	def test_overlap1(self):
		variant1 = Variant("chr1", 100, "AAAA", "A", 0)
		variant2 = Variant("chr1", 101, "A", "T", 1)
		variant3 = Variant("chr3", 100, "AAA", "A", 0)

		self.assertTrue(variants_overlap(variant1, variant2))
		self.assertTrue(variants_overlap(variant2, variant1))
		self.assertFalse(variants_overlap(variant1, variant3))
		self.assertFalse(variants_overlap(variant3, variant2))

	def test_overlap2(self):
		variant1 = Variant("chr1", 100, "AAAA", "A", 0)
		variant2 = Variant("chr1", 102, "AAAT", "A", 0)

		self.assertTrue(variants_overlap(variant1, variant2))
		self.assertTrue(variants_overlap(variant2, variant1))


	def test_overlap(self):
		variant1 = Variant("chr1", 100, "AA", "A", 0)
		variant2 = Variant("chr1", 102, "TT", "A", 0)

		self.assertFalse(variants_overlap(variant1, variant2))
		self.assertFalse(variants_overlap(variant2, variant1))


class Test_Cluster(unittest.TestCase):
	def test_cluster1(self):
		cluster = Cluster()
		variant1 = Variant("chr1", 100, "AAAA", "A", 0)
		variant2 = Variant("chr1", 101, "AAAA", "A", 0)

		cluster.add_var(variant1)
		cluster.add_var(variant2)

		variant3 = Variant("chr1", 103, "AAAA", "A", 0)
		self.assertTrue(cluster.overlaps_last(variant3))
		self.assertTrue(cluster.overlaps_last(variant2))
		self.assertTrue(cluster.overlaps_last(variant1))

	def test_cluster2(self):
		cluster = Cluster()
		variant1 = Variant("chr1", 100, "AAAA", "A", 0)

		self.assertFalse(cluster.overlaps_last(variant1))


class Test_group_variants(unittest.TestCase):
	def test_group1(self):
		variants = [ Variant("chr1", 0, "ATGCT", "A", 0),  Variant("chr1", 1, "T", "A", 1),  Variant("chr1", 3, "C", "A", 2),  Variant("chr1", 3, "CTAGC", "A", 3)]
		assignments = group_variants(variants)

		expected_assignments = [0, 1, 1, 2]
		self.assertEqual(assignments, expected_assignments)


	def test_group2(self):
		variants = [Variant("chr1", 0, "TAG", "A", 0), Variant("chr1", 2, "GCC", "A", 1), Variant("chr1", 4, "CTT", "A", 2)]
		assignments = group_variants(variants)

		expected_assignments = [0, 1, 0]
		self.assertEqual(assignments, expected_assignments)


	def test_group3(self):
		variants = [Variant("chr1", 0, "TAG", "A", 0), Variant("chr1", 1, "GCC", "A", 1), Variant("chr1", 2, "CTT", "A", 2)]
		assignments = group_variants(variants)

		expected_assignments = [0, 1, 2]
		self.assertEqual(assignments, expected_assignments)


class Test_run_paths(unittest.TestCase):
	def test_run_paths(self):

		vcf_lines = [
			"chr1	0	.	ATGCT	A	.	PASS	.	GT	1	0	0",
			"chr1	1	.	T	A	.	PASS	.	GT	0	1	0",
			"chr1	3	.	C	A	.	PASS	.	GT	0	1	0",
			"chr1	3	.	CTAGC	A	.	PASS	.	GT	0	0	1",
			"chr1	10	.	ATGATG	G	.	PASS	.	GT	1	0	0" 	]

		filename = "tmp.vcf.gz"
		with gzip.open(filename, 'wt') as tmp:
			for line in vcf_lines:
				tmp.write(line + '\n')
		run_paths(filename, None)


class Test_print_vcf(unittest.TestCase):
	def test_print_vcf(self):
		variants = [ Variant("chr1", 0, "ATGCT", "A", 0),  Variant("chr1", 1, "T", "A", 1),  Variant("chr1", 3, "C", "A", 2),  Variant("chr1", 3, "CTAGC", "A", 3)]
		assignments = group_variants(variants)
		print_vcf(variants, assignments)



if __name__ == '__main__':
	unittest.main()

