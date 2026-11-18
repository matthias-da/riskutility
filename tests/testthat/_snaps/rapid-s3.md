# print output is stable

    Code
      print(rapid(o, s, qi, "y", model_type = "lm", seed = 1))
    Output
      RAPID of 'y' given 2 quasi-identifiers (x1, x2), attacker lm: identity scale (chosen automatically)
      Rule: band: symmetric error < 10 %, gain > 0.3
      RAPID = 0.567: 85 of 150 records at risk
      Chance level = 0.279 (10 permutations)
      Excess = 0.288
      Next: summary(), plot(), confint(), rapid_test()

---

    Code
      print(rapid(o, s, qi, "y", model_type = "lm", num_eval_method = "point", seed = 1))
    Output
      RAPID of 'y' given 2 quasi-identifiers (x1, x2), attacker lm: identity scale (chosen automatically)
      Rule: point: symmetric error < 10 %
      RAPID = 0.767: 115 of 150 records at risk
      Chance level = 0.462 (10 permutations)
      Excess = 0.305
      Next: summary(), plot(), confint(), rapid_test()

---

    Code
      print(rapid(to_classes(o), to_classes(s, ref = o), qi, "y", model_type = "cart",
      seed = 1))
    Output
      RAPID of 'y' given 2 quasi-identifiers (x1, x2), attacker cart: 5 leaves, pruned by cross-validation
      Rule: RCS_marginal: gain over the marginal baseline > 0.3
      RAPID = 0.560: 84 of 150 records at risk
      Chance level = 0.341 (10 permutations)
      Excess = 0.219
      Next: summary(), plot(), confint(), rapid_test()

