# Code clean-up
- break up flowHist into smaller reusable chunks, to facilitate stepping through the init process, and redoing individual steps
- use consistent names for flowHist objects: pick one of fh, self etc and replace all the others.
- use consistent syntax - pass and return complete flowHist objects, rather than individual slots
- incorporate flowHist slots into flowFrame objects?

# Features
- linearity
- initial peak detection using FFT transform
- export results as a table for saving to disk

