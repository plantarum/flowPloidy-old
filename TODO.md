# Code clean-up
- [ ] break up flowHist into smaller reusable chunks, to facilitate stepping through the init process, and redoing individual steps
- [ ] use consistent names for flowHist objects: pick one of fh, self etc and replace all the others.
- [ ] use consistent syntax - pass and return complete flowHist objects, rather than individual slots
- [ ] incorporate flowHist slots into flowFrame objects?
- [x] use proper named lists for the components

# Features
- [ ] linearity
- [ ] initial peak detection with smoothing using kza::kz (or FFT?)
- [x] export results as a table for saving to disk
- [ ] save and restore parameter estimates from nls
- [x] batch processing
- [ ] switch to nlsLM for regressions


