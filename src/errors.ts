export class BelowIntrinsicError extends Error {
  constructor(message = "The option price is below intrinsic value.") {
    super(message);
    this.name = "BelowIntrinsicError";
  }
}

export class AboveMaximumError extends Error {
  constructor(message = "The option price is above the maximum option value.") {
    super(message);
    this.name = "AboveMaximumError";
  }
}
