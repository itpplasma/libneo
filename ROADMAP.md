# libneo roadmap

Status: demand-driven modernization plan, 2026-10-05.

libneo is mature shared infrastructure with historical layers. The objective is
not a clean-slate rewrite. The objective is to make every newly exercised path
more correct, testable and reusable while preserving the accumulated
format/convention knowledge on which ITP plasma codes depend.

## L0 — preserve shared production value

Keep current NEO-2, SIMPLE, MEPHIT, TIAGO and other downstream use cases
working. Correctness fixes take priority over source compatibility when the old
behavior is physically wrong, but provide explicit migration notes or bounded
compatibility shims when practical.

Every convention-sensitive reader/converter should accumulate independent
analytic, round-trip or cross-code oracles rather than relying only on another
consumer producing plausible output.

## L1 — KIN6D first-use hardening gate

As soon as a KIN6D milestone depends on a libneo capability, audit that exact
path before treating it as production infrastructure.

For the touched reader/converter/evaluator:

1. capture the KIN6D case as a minimal libneo regression where possible;
2. verify units, signs, field periods, flux normalization, orientation and
   coordinate conventions;
3. fix interpolation/accuracy/caching defects upstream;
4. remove or encapsulate global state when it blocks reentrancy, threading,
   multiple simultaneous instances or derivative testing;
5. benchmark only if the path is material to the KIN6D workload;
6. re-run relevant existing downstream tests.

Do not create a KIN6D-private corrected GEQDSK, VMEC, Boozer, coil or field
implementation when the correction is generic.

## L2 — optional differentiable field capabilities

Add stronger interfaces only when concrete consumers appear. Keep the simple
value interface usable.

Candidate reusable capabilities include:

- spatial field Jacobians and directional derivatives;
- consistent derivatives of coordinate/interpolation maps;
- batch JVP-like evaluation where it improves a real consumer;
- derivative regression tests against analytic/complex-step/finite-difference
  oracles.

Do not put physical equilibrium parameter adjoints here. Differentiating
(R(U,\theta)=0) and its implicit solve belongs to KIN6D (or another owning
physics solver), not the import library.

## L3 — representation-level error information

For imported/interpolated data, add error information where it can be stated
generically and checked:

- interpolation residual/error indicators;
- spectral/truncation metadata;
- consistency residuals such as divergence/Ampere/format identities;
- interval/ball evaluation for a representation when a concrete certifier
  needs it.

These bounds describe libneo's representation/evaluation, not the error of an
MHD/DK/GK/FK physical model.

If retrofitting rigorous evaluation into a historical backend is
disproportionate, normalize/export its data into a canonical representation
that KIN6D or another consumer can differentiate and certify independently.

## L4 — reduce historical coupling only on contact

When a touched path has hidden global state, unnecessary transitive
dependencies, compiler coupling or duplicated generic utilities, refactor only
as far as the current consumers justify. Prefer explicit state and small
capability interfaces.

Generic functionality that clearly belongs to FortIO, FortNum, FortFEM or
another controlled library may migrate there after at least two real consumers
establish the stable abstraction. Do not reorganize the repository merely to
make the directory tree look modern.

## Permanent ecosystem rule

A KIN6D-driven libneo improvement is an ecosystem improvement.

The normal flow is:

    KIN6D/TIAGO/SIMPLE/MEPHIT/... exposes need
        -> reproduce in libneo
        -> test and repair/refactor libneo
        -> verify relevant downstream consumer
        -> update/pin consumers

The downstream adapter remains narrow and replaceable, but production
mathematics is not duplicated simply to preserve that replaceability.

Chris&AI
