import { erfcxCody, inverseNormCdf, normCdf, normPdf } from "cody-special";
import {
  convexRationalCubicControlParameterToFitSecondDerivativeAtLeftSide,
  convexRationalCubicControlParameterToFitSecondDerivativeAtRightSide,
  rationalCubicInterpolation
} from "piecewise-rational";
import { AboveMaximumError, BelowIntrinsicError } from "./errors.js";
import {
  DBL_EPSILON,
  DBL_MAX,
  DBL_MIN,
  DENORMALIZATION_CUTOFF,
  FOURTH_ROOT_DBL_EPSILON,
  ONE_OVER_SQRT_TWO,
  ONE_OVER_SQRT_TWO_PI,
  PI_OVER_SIX,
  SIXTEENTH_ROOT_DBL_EPSILON,
  SQRT_DBL_MAX,
  SQRT_DBL_MIN,
  SQRT_ONE_OVER_THREE,
  SQRT_PI_OVER_TWO,
  SQRT_THREE,
  SQRT_TWO_PI,
  TWO_PI,
  TWO_PI_OVER_SQRT_TWENTY_SEVEN
} from "./constants.js";

export { normCdf, normPdf } from "cody-special";
export { erfCody, erfcCody, erfcxCody, inverseNormCdf } from "cody-special";
export * from "piecewise-rational";
export * from "./errors.js";

export const impliedVolatilityMaximumIterations = 2;
const asymptoticExpansionAccuracyThreshold = -10;
const smallTExpansionOfNormalizedBlackThreshold = 2 * SIXTEENTH_ROOT_DBL_EPSILON;

function square(x: number): number {
  return x * x;
}

function isBelowHorizon(x: number): boolean {
  return Math.abs(x) < DENORMALIZATION_CUTOFF;
}

function householderFactor(newton: number, halley: number, hh3: number): number {
  return (1 + 0.5 * halley * newton) / (1 + newton * (halley + hh3 * newton / 6));
}

function computeFLowerMapAndFirstTwoDerivatives(x: number, s: number): [number, number, number] {
  const ax = Math.abs(x);
  const z = SQRT_ONE_OVER_THREE * ax / s;
  const y = z * z;
  const s2 = s * s;
  const Phi = normCdf(-z);
  const phi = normPdf(z);
  const fpp = PI_OVER_SIX * y / (s2 * s) * Phi *
    (8 * SQRT_THREE * s * ax + (3 * s2 * (s2 - 8) - 8 * x * x) * Phi / phi) *
    Math.exp(2 * y + 0.25 * s2);
  if (isBelowHorizon(s)) return [0, 1, fpp];
  const Phi2 = Phi * Phi;
  const fp = TWO_PI * y * Phi2 * Math.exp(y + 0.125 * s * s);
  const f = isBelowHorizon(x) ? 0 : TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * (Phi2 * Phi);
  return [f, fp, fpp];
}

function computeFUpperMapAndFirstTwoDerivatives(x: number, s: number): [number, number, number] {
  const f = normCdf(-0.5 * s);
  if (isBelowHorizon(x)) return [f, -0.5, 0];
  const w = square(x / s);
  const fp = -0.5 * Math.exp(0.5 * w);
  const fpp = SQRT_PI_OVER_TWO * Math.exp(w + 0.125 * s * s) * w / s;
  return [f, fp, fpp];
}

function inverseFLowerMap(x: number, f: number): number {
  return isBelowHorizon(f) ? 0 : Math.abs(x / (SQRT_THREE * inverseNormCdf(Math.pow(f / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * Math.abs(x)), 1 / 3))));
}

function inverseFUpperMap(f: number): number {
  return -2 * inverseNormCdf(f);
}

function normalizedIntrinsic(x: number, q: number): number {
  if (q * x <= 0) return 0;
  const sign = q < 0 ? -1 : 1;
  const x2 = x * x;
  if (x2 < 98 * FOURTH_ROOT_DBL_EPSILON) {
    return Math.abs(Math.max(sign * x * (1 + x2 * (1 / 24 + x2 * (1 / 1920 + x2 * (1 / 322560 + x2 / 92897280)))), 0));
  }
  const bMax = Math.exp(0.5 * x);
  return Math.abs(Math.max(sign * (bMax - 1 / bMax), 0));
}

export function normalisedIntrinsic(x: number, q: number): number {
  return normalizedIntrinsic(x, q);
}

export function normalizedIntrinsicCall(x: number): number {
  return normalizedIntrinsic(x, 1);
}

export const normalisedIntrinsicCall = normalizedIntrinsicCall;

function normalizedBlackCallUsingNormCdf(x: number, s: number): number {
  if (s <= 0) return normalizedIntrinsicCall(x);
  const h = x / s;
  const t = 0.5 * s;
  const bMax = Math.exp(0.5 * x);
  const b = normCdf(h + t) * bMax - normCdf(h - t) / bMax;
  return Math.abs(Math.max(b, 0));
}

function smallTExpansionOfNormalizedBlackCall(h: number, t: number): number {
  const a = 1 + h * (0.5 * SQRT_TWO_PI) * erfcxCody(-ONE_OVER_SQRT_TWO * h);
  const w = t * t;
  const h2 = h * h;
  const c1 = (-1 + 3 * a + a * h2) / 6;
  const c2 = (-7 + 15 * a + h2 * (-1 + 10 * a + a * h2)) / 120;
  const c3 = (-57 + 105 * a + h2 * (-18 + 105 * a + h2 * (-1 + 21 * a + a * h2))) / 5040;
  const c4 = (-561 + 945 * a + h2 * (-285 + 1260 * a + h2 * (-33 + 378 * a + h2 * (-1 + 36 * a + a * h2)))) / 362880;
  const c5 = (-6555 + 10395 * a + h2 * (-4680 + 17325 * a + h2 * (-840 + 6930 * a + h2 * (-52 + 990 * a + h2 * (-1 + 55 * a + a * h2))))) / 39916800;
  const c6 = (-89055 + 135135 * a + h2 * (-82845 + 270270 * a + h2 * (-20370 + 135135 * a + h2 * (-1926 + 25740 * a + h2 * (-75 + 2145 * a + h2 * (-1 + 78 * a + a * h2)))))) / 6227020800;
  const expansion = 2 * t * (a + w * (c1 + w * (c2 + w * (c3 + w * (c4 + w * (c5 + c6 * w))))));
  const b = ONE_OVER_SQRT_TWO_PI * Math.exp(-0.5 * (h * h + t * t)) * expansion;
  return Math.abs(Math.max(b, 0));
}

function normalizedBlackCallUsingErfcx(h: number, t: number): number {
  const b = 0.5 * Math.exp(-0.5 * (h * h + t * t)) *
    (erfcxCody(-ONE_OVER_SQRT_TWO * (h + t)) - erfcxCody(-ONE_OVER_SQRT_TWO * (h - t)));
  return Math.abs(Math.max(b, 0));
}

export function normalizedBlackCall(x: number, s: number): number {
  if (x > 0) return normalizedIntrinsicCall(x) + normalizedBlackCall(-x, s);
  const ax = Math.abs(x);
  if (s <= ax * DENORMALIZATION_CUTOFF) return normalizedIntrinsicCall(x);
  if (
    x < s * asymptoticExpansionAccuracyThreshold &&
    0.5 * s * s + x < s * (smallTExpansionOfNormalizedBlackThreshold + asymptoticExpansionAccuracyThreshold)
  ) {
    return normalizedBlackCallUsingErfcx(x / s, 0.5 * s);
  }
  if (0.5 * s < smallTExpansionOfNormalizedBlackThreshold) {
    return smallTExpansionOfNormalizedBlackCall(x / s, 0.5 * s);
  }
  if (x + 0.5 * s * s > s * 0.85) {
    return normalizedBlackCallUsingNormCdf(x, s);
  }
  return normalizedBlackCallUsingErfcx(x / s, 0.5 * s);
}

export const normalisedBlackCall = normalizedBlackCall;

export function normalizedBlack(x: number, s: number, q: number): number {
  return normalizedBlackCall(q < 0 ? -x : x, s);
}

export const normalisedBlack = normalizedBlack;

export function black(F: number, K: number, sigma: number, T: number, q: number): number {
  const intrinsic = Math.abs(Math.max(q < 0 ? K - F : F - K, 0));
  if (q * (F - K) > 0) return intrinsic + black(F, K, sigma, T, -q);
  return Math.max(intrinsic, Math.sqrt(F) * Math.sqrt(K) * normalizedBlack(Math.log(F / K), sigma * Math.sqrt(T), q));
}

export function normalizedVega(x: number, s: number): number {
  const ax = Math.abs(x);
  if (ax <= 0) return ONE_OVER_SQRT_TWO_PI * Math.exp(-0.125 * s * s);
  return s <= 0 || s <= ax * SQRT_DBL_MIN ? 0 : ONE_OVER_SQRT_TWO_PI * Math.exp(-0.5 * (square(x / s) + square(0.5 * s)));
}

export const normalisedVega = normalizedVega;

function uncheckedNormalizedImpliedVolatility(beta: number, x: number, q: number, iterations: number): number {
  if (q * x > 0) {
    beta = Math.abs(Math.max(beta - normalizedIntrinsic(x, q), 0));
    q = -q;
  }
  if (q < 0) {
    x = -x;
    q = -q;
  }
  if (beta <= 0 || beta < DENORMALIZATION_CUTOFF) return 0;
  const bMax = Math.exp(0.5 * x);
  if (beta >= bMax) throw new AboveMaximumError();

  let count = 0;
  let directionReversalCount = 0;
  let f = -DBL_MAX;
  let s = -DBL_MAX;
  let ds = s;
  let dsPrevious = 0;
  let sLeft = DBL_MIN;
  let sRight = DBL_MAX;
  const sC = Math.sqrt(Math.abs(2 * x));
  const bC = normalizedBlackCall(x, sC);
  const vC = normalizedVega(x, sC);

  if (beta < bC) {
    const sL = sC - bC / vC;
    const bL = normalizedBlackCall(x, sL);
    if (beta < bL) {
      const [fLowerMapL, dFLowerMapLDBeta, d2FLowerMapLDBeta2] = computeFLowerMapAndFirstTwoDerivatives(x, sL);
      const rLL = convexRationalCubicControlParameterToFitSecondDerivativeAtRightSide(0, bL, 0, fLowerMapL, 1, dFLowerMapLDBeta, d2FLowerMapLDBeta2, true);
      f = rationalCubicInterpolation(beta, 0, bL, 0, fLowerMapL, 1, dFLowerMapLDBeta, rLL);
      if (!(f > 0)) {
        const t = beta / bL;
        f = (fLowerMapL * t + bL * (1 - t)) * t;
      }
      s = inverseFLowerMap(x, f);
      sRight = sL;
      while (count < iterations && Math.abs(ds) > DBL_EPSILON * s) {
        if (ds * dsPrevious < 0) directionReversalCount += 1;
        if (count > 0 && (directionReversalCount === 3 || !(s > sLeft && s < sRight))) {
          s = 0.5 * (sLeft + sRight);
          if (sRight - sLeft <= DBL_EPSILON * s) break;
          directionReversalCount = 0;
          ds = 0;
        }
        dsPrevious = ds;
        const b = normalizedBlackCall(x, s);
        const bp = normalizedVega(x, s);
        if (b > beta && s < sRight) sRight = s;
        else if (b < beta && s > sLeft) sLeft = s;
        if (b <= 0 || bp <= 0) {
          ds = 0.5 * (sLeft + sRight) - s;
        } else {
          const lnB = Math.log(b);
          const lnBeta = Math.log(beta);
          const bpob = bp / b;
          const h = x / s;
          const bHalley = h * h / s - s / 4;
          const newton = (lnBeta - lnB) * lnB / lnBeta / bpob;
          const halley = bHalley - bpob * (1 + 2 / lnB);
          const bHh3 = bHalley * bHalley - 3 * square(h / s) - 0.25;
          const hh3 = bHh3 + 2 * square(bpob) * (1 + 3 / lnB * (1 + 1 / lnB)) - 3 * bHalley * bpob * (1 + 2 / lnB);
          ds = newton * householderFactor(newton, halley, hh3);
        }
        ds = Math.max(-0.5 * s, ds);
        s += ds;
        count += 1;
      }
      return s;
    }
    const vL = normalizedVega(x, sL);
    const rLM = convexRationalCubicControlParameterToFitSecondDerivativeAtRightSide(bL, bC, sL, sC, 1 / vL, 1 / vC, 0, false);
    s = rationalCubicInterpolation(beta, bL, bC, sL, sC, 1 / vL, 1 / vC, rLM);
    sLeft = sL;
    sRight = sC;
  } else {
    const sH = vC > DBL_MIN ? sC + (bMax - bC) / vC : sC;
    const bH = normalizedBlackCall(x, sH);
    if (beta <= bH) {
      const vH = normalizedVega(x, sH);
      const rHM = convexRationalCubicControlParameterToFitSecondDerivativeAtLeftSide(bC, bH, sC, sH, 1 / vC, 1 / vH, 0, false);
      s = rationalCubicInterpolation(beta, bC, bH, sC, sH, 1 / vC, 1 / vH, rHM);
      sLeft = sC;
      sRight = sH;
    } else {
      const [fUpperMapH, dFUpperMapHDBeta, d2FUpperMapHDBeta2] = computeFUpperMapAndFirstTwoDerivatives(x, sH);
      if (d2FUpperMapHDBeta2 > -SQRT_DBL_MAX && d2FUpperMapHDBeta2 < SQRT_DBL_MAX) {
        const rHH = convexRationalCubicControlParameterToFitSecondDerivativeAtLeftSide(bH, bMax, fUpperMapH, 0, dFUpperMapHDBeta, -0.5, d2FUpperMapHDBeta2, true);
        f = rationalCubicInterpolation(beta, bH, bMax, fUpperMapH, 0, dFUpperMapHDBeta, -0.5, rHH);
      }
      if (f <= 0) {
        const h = bMax - bH;
        const t = (beta - bH) / h;
        f = (fUpperMapH * (1 - t) + 0.5 * h * t) * (1 - t);
      }
      s = inverseFUpperMap(f);
      sLeft = sH;
      if (beta > 0.5 * bMax) {
        while (count < iterations && Math.abs(ds) > DBL_EPSILON * s) {
          if (ds * dsPrevious < 0) directionReversalCount += 1;
          if (count > 0 && (directionReversalCount === 3 || !(s > sLeft && s < sRight))) s = 0.5 * (sLeft + sRight);
          if (sRight - sLeft <= DBL_EPSILON * s) break;
          directionReversalCount = 0;
          ds = 0;
          dsPrevious = ds;
          const b = normalizedBlackCall(x, s);
          const bp = normalizedVega(x, s);
          if (b > beta && s < sRight) sRight = s;
          else if (b < beta && s > sLeft) sLeft = s;
          if (b >= bMax || bp <= DBL_MIN) {
            ds = 0.5 * (sLeft + sRight) - s;
          } else {
            const bMaxMinusB = bMax - b;
            const g = Math.log((bMax - beta) / bMaxMinusB);
            const gp = bp / bMaxMinusB;
            const bHalley = square(x / s) / s - s / 4;
            const bHh3 = bHalley * bHalley - 3 * square(x / (s * s)) - 0.25;
            const newton = -g / gp;
            const halley = bHalley + gp;
            const hh3 = bHh3 + gp * (2 * gp + 3 * bHalley);
            ds = newton * householderFactor(newton, halley, hh3);
          }
          ds = Math.max(-0.5 * s, ds);
          s += ds;
          count += 1;
        }
        return s;
      }
    }
  }

  while (count < iterations && Math.abs(ds) > DBL_EPSILON * s) {
    if (ds * dsPrevious < 0) directionReversalCount += 1;
    if (count > 0 && (directionReversalCount === 3 || !(s > sLeft && s < sRight))) {
      s = 0.5 * (sLeft + sRight);
      if (sRight - sLeft <= DBL_EPSILON * s) break;
      directionReversalCount = 0;
      ds = 0;
    }
    dsPrevious = ds;
    const b = normalizedBlackCall(x, s);
    const bp = normalizedVega(x, s);
    if (b > beta && s < sRight) sRight = s;
    else if (b < beta && s > sLeft) sLeft = s;
    const newton = (beta - b) / bp;
    const halley = square(x / s) / s - s / 4;
    const hh3 = halley * halley - 3 * square(x / (s * s)) - 0.25;
    ds = Math.max(-0.5 * s, newton * householderFactor(newton, halley, hh3));
    s += ds;
    count += 1;
  }
  return s;
}

export function normalizedImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations(beta: number, x: number, q: number, iterations: number): number {
  if (q * x > 0) {
    beta -= normalizedIntrinsic(x, q);
    q = -q;
  }
  if (beta < 0) throw new BelowIntrinsicError();
  return uncheckedNormalizedImpliedVolatility(beta, x, q, iterations);
}

export const normalisedImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations =
  normalizedImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations;

export function normalizedImpliedVolatilityFromATransformedRationalGuess(beta: number, x: number, q: number): number {
  return normalizedImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations(beta, x, q, impliedVolatilityMaximumIterations);
}

export const normalisedImpliedVolatilityFromATransformedRationalGuess = normalizedImpliedVolatilityFromATransformedRationalGuess;

export function impliedVolatilityFromATransformedRationalGuessWithLimitedIterations(price: number, F: number, K: number, T: number, q: number, iterations: number): number {
  const intrinsic = Math.abs(Math.max(q < 0 ? K - F : F - K, 0));
  if (price < intrinsic) throw new BelowIntrinsicError();
  const maxPrice = q < 0 ? K : F;
  if (price >= maxPrice) throw new AboveMaximumError();
  let x = Math.log(F / K);
  if (q * x > 0) {
    price = Math.abs(Math.max(price - intrinsic, 0));
    q = -q;
  }
  return uncheckedNormalizedImpliedVolatility(price / (Math.sqrt(F) * Math.sqrt(K)), x, q, iterations) / Math.sqrt(T);
}

export function impliedVolatilityFromATransformedRationalGuess(price: number, F: number, K: number, T: number, q: number): number {
  return impliedVolatilityFromATransformedRationalGuessWithLimitedIterations(price, F, K, T, q, impliedVolatilityMaximumIterations);
}
