import unittest

from test_structure import TestStructure


def build_suite():
    # Master
    master_suite = unittest.TestSuite()

    # Structure
    suite1 = unittest.TestSuite()
    suite1.addTest(unittest.makeSuite(TestStructure))

    # Add tests
    master_suite.addTest(suite1)

    return master_suite


if __name__ == '__main__':
    suite = build_suite()
    unittest.TextTestRunner(verbosity=2).run(suite)
