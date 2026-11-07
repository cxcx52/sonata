"""Exercise a valid slope, a short-kernel fork, and rejected input records."""
from dataclasses import replace
from random import Random
import unittest

from extraction.ball_challenges import BallIndex, FilteredBall, filtered_fiber_ranks, negacyclic_mul
from extraction.coordinate import DEGREE, Q, Record, compare_coordinate


def poly(*coefficients):
    return tuple(coefficients) + (0,) * (DEGREE - len(coefficients))


class CoordinateCheck(unittest.TestCase):
    def setUp(self):
        self.matrix = ((poly(1), poly(1)),)
        self.slope = (poly(3, 1), poly(-1, 2))
        self.image = (poly(2, 3),)
        self.offset = (poly(0),)
        self.records = tuple(Record(c, tuple(negacyclic_mul(c,w) for w in self.slope))
                             for c in (poly(1), poly(0,1), poly(-1,1)))

    def evaluate(self, records):
        return compare_coordinate(self.matrix,self.image,self.offset,records,100)

    def test_slope(self):
        got = self.evaluate(self.records)
        self.assertEqual(got.kind, "slope")
        self.assertEqual(got.vector,tuple(tuple(y % Q for y in x) for x in self.slope))

    def test_kernel(self):
        first, second, fork = self.records
        left, right = fork.response
        changed = ((left[0]+1,*left[1:]), (right[0]-1,*right[1:]))
        got = self.evaluate((first, second, replace(fork,response=changed)))
        self.assertEqual(got.kind,"kernel")
        self.assertTrue(any(v for x in got.vector for v in x))
        self.assertTrue(all((l+r) % Q == 0 for l,r in zip(*got.vector)))
        self.assertLessEqual(sum(v*v for x in got.vector for v in x), (8*26)**2*100)

    def test_rejections(self):
        first, second, fork = self.records
        for bad in (first, replace(fork,response=(poly(0),poly(0))), replace(fork,challenge=poly(12))):
            with self.assertRaises(ValueError):
                self.evaluate((first,second,bad))

    def test_pivot_fork_support(self):
        ball = BallIndex(4, 2)
        domain = FilteredBall(ball, 2)
        center = ball.rank((1,0,0,0))
        ranks = list(filtered_fiber_ranks(domain, center, Random(1).randrange, nonzero=True))
        expected = {i for i in range(ball.size) if i not in (center,ball.zero_rank)
                    and domain.member(ball.unrank(i))}
        self.assertEqual(set(ranks), expected)
        self.assertEqual(len(ranks), len(expected))


if __name__ == "__main__":
    unittest.main()
