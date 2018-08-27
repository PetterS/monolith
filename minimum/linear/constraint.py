class Constraint:
	def __init__(self, lb, sum, ub) -> None:
		"""Creates the constraint lb ≤ sum ≤ ub."""
		self.lb = lb
		self.sum = sum.copy()
		self.ub = ub

		if self.lb is not None:
			self.lb -= self.sum.constant
		if self.ub is not None:
			self.ub -= self.sum.constant
		self.sum.constant = 0

		if len(self.sum.coefficients) == 0:
			assert self.lb is None or self.lb <= 0, "Can not add constraint impossible to fulfill."
			assert self.ub is None or 0 <= self.ub, "Can not add constraint impossible to fulfill."

	def __str__(self) -> str:
		result = ""
		if self.lb is not None:
			result += str(self.lb) + " ≤ "
		result += str(self.sum)
		if self.ub is not None:
			result += " ≤ " + str(self.ub)
		return result
