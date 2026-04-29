import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import {
  AboveMaximumError,
  BelowIntrinsicError,
  black,
  impliedVolatilityFromATransformedRationalGuess,
  impliedVolatilityFromATransformedRationalGuessWithLimitedIterations,
  normCdf,
  normalisedBlack,
  normalisedBlackCall,
  normalisedImpliedVolatilityFromATransformedRationalGuess,
  normalisedImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations,
  normalisedVega
} from "../src/index.js";

function close(actual: number, expected: number, tolerance = 1e-10): void {
  assert.ok(Math.abs(actual - expected) <= tolerance || Math.abs(actual - expected) <= tolerance * Math.abs(expected), `${actual} != ${expected}`);
}

test("normalised Black and Black price match Python reference samples", () => {
  close(normCdf(0), 0.5, 1e-15);
  close(normalisedBlackCall(0, 0.2), 0.07965567455405798, 1e-14);
  close(normalisedBlack(0, 0.2, 1), 0.07965567455405798, 1e-14);
  close(black(100, 100, 0.2, 1, 1), 7.965567455405798, 1e-12);
  close(normalisedVega(0, 0.2), 0.3969525474770118, 1e-14);
});

test("implied volatility inverts normalized and standard Black prices", () => {
  const beta = normalisedBlack(0.1, 0.25, 1);
  close(normalisedImpliedVolatilityFromATransformedRationalGuess(beta, 0.1, 1), 0.25, 1e-12);
  const price = black(105, 100, 0.3, 2, 1);
  close(impliedVolatilityFromATransformedRationalGuess(price, 105, 100, 2, 1), 0.3, 1e-12);
});

test("implied volatility rejects prices outside bounds", () => {
  assert.throws(() => impliedVolatilityFromATransformedRationalGuess(9.99, 110, 100, 1, 1), BelowIntrinsicError);
  assert.throws(() => impliedVolatilityFromATransformedRationalGuess(100, 110, 100, 1, -1), AboveMaximumError);
});

const regressionCases = [
  {
    input: { q: 1, s: 0.3384896249542857, T: 0.89, F: 111, x: 0.5677720320717676, K: 246, z: 0.626639629966119, sigma: 0.6907245366960202, N: 3 },
    output: { normCdf: 0.7345522576674526, black: 5.578476100968555, normalisedBlack: 0.5818992767215997, normalisedBlackCall: 0.5818992767215997, normalisedVega: 0.09632340394779765, normalisedIv: 0.33848962495428553, normalisedIvLimited: 0.3384896249542856, iv: 0.6907245366960202, ivLimited: 0.6907245366960202 }
  },
  {
    input: { q: 1, s: 0.5170014103772701, T: 0.71, F: 442, x: 0.6527949377834605, K: 481, z: 0.44154194429522964, sigma: 0.656371761477067, N: 0 },
    output: { normCdf: 0.6705896485551225, black: 82.1858144664737, normalisedBlack: 0.6893418675139219, normalisedBlackCall: 0.6893418675139219, normalisedVega: 0.17386086261690376, normalisedIv: 0.5170014103772699, normalisedIvLimited: 0.5170219815983038, iv: 0.6563717614770669, ivLimited: 0.656255585693245 }
  },
  {
    input: { q: -1, s: 0.20945347091101296, T: 1.35, F: 340, x: 0.4130360224654511, K: 114, z: 0.8882590309465282, sigma: 0.7797391097354639, N: 1 },
    output: { normCdf: 0.8127992840462526, black: 9.222103601508342, normalisedBlack: 0.0019082840985949082, normalisedBlackCall: 0.4178865501642314, normalisedVega: 0.05676951548164378, normalisedIv: 0.20945347091101296, normalisedIvLimited: 0.20945349849184025, iv: 0.7797391097354639, ivLimited: 0.7797391097354656 }
  }
] as const;

test("selected LetsBeRational regression grid values", () => {
  for (const { input, output } of regressionCases) {
    const beta = normalisedBlack(input.x, input.s, input.q);
    close(normCdf(input.z), output.normCdf, 1e-12);
    close(black(input.F, input.K, input.sigma, input.T, input.q), output.black, 1e-10);
    close(beta, output.normalisedBlack, 1e-12);
    close(normalisedBlackCall(input.x, input.s), output.normalisedBlackCall, 1e-12);
    close(normalisedVega(input.x, input.s), output.normalisedVega, 1e-12);
    close(normalisedImpliedVolatilityFromATransformedRationalGuess(beta, input.x, input.q), output.normalisedIv, 1e-10);
    close(
      normalisedImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations(beta, input.x, input.q, input.N),
      output.normalisedIvLimited,
      1e-10
    );
    const price = black(input.F, input.K, input.sigma, input.T, input.q);
    close(impliedVolatilityFromATransformedRationalGuess(price, input.F, input.K, input.T, input.q), output.iv, 1e-10);
    close(
      impliedVolatilityFromATransformedRationalGuessWithLimitedIterations(price, input.F, input.K, input.T, input.q, input.N),
      output.ivLimited,
      1e-10
    );
  }
});

type GridInput = {
  q: number;
  s: number;
  T: number;
  F: number;
  x: number;
  K: number;
  z: number;
  sigma: number;
  N: number;
};

type GridOutput = {
  norm_cdf: number;
  black: number;
  normalised_black: number;
  normalised_black_call: number;
  normalised_vega: number;
  normalised_implied_volatility_from_a_transformed_rational_guess: number;
  normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations: number;
  implied_volatility_from_a_transformed_rational_guess: number | null;
  implied_volatility_from_a_transformed_rational_guess_with_limited_iterations: number | null;
};

const grid = JSON.parse(
  readFileSync("test/TestValues.json", "utf8").replace(/\bNaN\b/g, "null")
) as { input: GridInput[]; output: GridOutput[] };

type CppVectors = {
  tolerancePolicy: {
    absolute: number;
    relative: number;
  };
  cases: Array<{
    id: string;
    input: {
      x: number;
      s?: number;
      theta: number;
      priceFractionOfBMax?: number;
    };
    output: {
      bmax: number;
      normalisedBlack?: number;
      normalisedVega?: number;
      normalisedImpliedVolatilityFromPrice?: number;
      normalisedPrice?: number;
      normalisedImpliedVolatility?: number;
    };
  }>;
};

const cppVectors = JSON.parse(
  readFileSync("node_modules/vollib-test-vectors/vectors/lets-be-rational/binary64/normalised-black-reduced-domain.json", "utf8")
) as CppVectors;

test("full LetsBeRational regression value grid", () => {
  for (let i = 0; i < grid.input.length; i += 1) {
    const input = grid.input[i];
    const output = grid.output[i];
    const beta = normalisedBlack(input.x, input.s, input.q);

    close(normCdf(input.z), output.norm_cdf, 1e-12);
    close(black(input.F, input.K, input.sigma, input.T, input.q), output.black, 1e-10);
    close(beta, output.normalised_black, 1e-12);
    close(normalisedBlackCall(input.x, input.s), output.normalised_black_call, 1e-12);
    close(normalisedVega(input.x, input.s), output.normalised_vega, 1e-12);
    close(
      normalisedImpliedVolatilityFromATransformedRationalGuess(beta, input.x, input.q),
      output.normalised_implied_volatility_from_a_transformed_rational_guess,
      1e-10
    );
    close(
      normalisedImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations(beta, input.x, input.q, input.N),
      output.normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations,
      1e-10
    );

    if (input.T > 0 && output.implied_volatility_from_a_transformed_rational_guess !== null) {
      const price = black(input.F, input.K, input.sigma, input.T, input.q);
      close(
        impliedVolatilityFromATransformedRationalGuess(price, input.F, input.K, input.T, input.q),
        output.implied_volatility_from_a_transformed_rational_guess,
        1e-9
      );
      close(
        impliedVolatilityFromATransformedRationalGuessWithLimitedIterations(price, input.F, input.K, input.T, input.q, input.N),
        output.implied_volatility_from_a_transformed_rational_guess_with_limited_iterations!,
        1e-9
      );
    }
  }
});

test("normalised Black matches shared C++ binary64 reduced-domain vectors", () => {
  const directTolerance = Math.max(cppVectors.tolerancePolicy.absolute, cppVectors.tolerancePolicy.relative);
  const impliedVolatilityTolerance = 2e-3;
  for (const row of cppVectors.cases) {
    const { x, theta } = row.input;
    assert.equal(theta, 1, `${row.id}: expected reduced-domain theta = +1`);
    assert.ok(x <= 0, `${row.id}: expected reduced-domain x <= 0`);

    if (row.input.s !== undefined) {
      const beta = normalisedBlack(x, row.input.s, theta);
      close(beta, row.output.normalisedBlack!, directTolerance);
      close(normalisedVega(x, row.input.s), row.output.normalisedVega!, directTolerance);

      if (row.output.normalisedImpliedVolatilityFromPrice! > 0 && beta < row.output.bmax) {
        close(
          normalisedImpliedVolatilityFromATransformedRationalGuess(beta, x, theta),
          row.output.normalisedImpliedVolatilityFromPrice!,
          impliedVolatilityTolerance
        );
      }
    } else {
      const beta = row.output.normalisedPrice!;
      close(
        normalisedImpliedVolatilityFromATransformedRationalGuess(beta, x, theta),
        row.output.normalisedImpliedVolatility!,
        impliedVolatilityTolerance
      );
    }
  }
});

test("normalised implied volatility round trips back to the shared C++ price", () => {
  const priceRoundTripTolerance = Math.max(cppVectors.tolerancePolicy.absolute, cppVectors.tolerancePolicy.relative);
  for (const row of cppVectors.cases) {
    const { x, theta } = row.input;
    const beta = row.input.s === undefined ? row.output.normalisedPrice! : row.output.normalisedBlack!;

    if (beta <= 0 || beta >= row.output.bmax) continue;

    const impliedVolatility = normalisedImpliedVolatilityFromATransformedRationalGuess(beta, x, theta);
    const roundTrippedPrice = normalisedBlack(x, impliedVolatility, theta);
    close(roundTrippedPrice, beta, priceRoundTripTolerance);
  }
});
