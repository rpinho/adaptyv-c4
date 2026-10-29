# Invalidated diagnostic run — do not select these variants

This first smoke run used the original installer protocol: bound input locally
repacked at mutation sites, but InterfaceAnalyzer's separated partners globally
repacked. The reference dG was **+626.91 REU**, so the relative ddG filter alone
could call variants retained while absolute binding energy was grossly unfavorable.
These outputs are retained to document the error; **none is a validated binding
candidate**. The current campaign uses shared whole-complex side-chain preparation
and also requires negative reference and variant dG before retaining a twin.

See `../energy_diagnostic.log`, `../prepack_diagnostic.log`, `../VALIDATION.md`
and the replacement `../campaign/` results. The unprepared source had dG +804.75;
input packing alone changed it to +8.93; common full repacking plus fixed-backbone
side-chain minimisation gave -44.66. This demonstrates sensitivity to preparation,
not experimental binding, and requires matched references throughout.
