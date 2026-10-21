# argon2id-spark: Argon2id in SPARK, on a proved BLAKE2b

*Project document, 21 October 2026. "argon2id-spark" is a working name
(crate `argon2id_spark`), to be checked before publication.*

> **Status of every claim in this document.**
> - **Read at the source**: RFC 9106 (Argon2), sections cited as §n of the
>   RFC.
> - **Measured** (21 October 2026, Linux x86-64): the `blake2b` crate in
>   the Alire index (0.1.1, BSD-3-Clause, "in SPARK, proven against a tested
>   specification"); no Argon2 crate in the index; the reference
>   implementation installed (`argon2` command line and `libargon2`
>   20190702 with its header); gnatprove FSF 16.1.0.
> - **[to verify]**: the `blake2b` crate's API (variable output length,
>   streaming), the PHC string format's specification, GNAT's support for
>   Ada 2022 parallel constructs, and anything about vector instructions
>   from Ada.

---

## 1. Summary

Argon2id is the password-hashing and key-derivation function RFC 9106
recommends. Ada has no implementation (measured: no Argon2 crate in the
Alire index); an Ada program binds libargon2 in C, or uses OpenSSL's KDF
where its version has one.

argon2id-spark is Argon2id (and Argon2i and Argon2d, which share the code)
written in **SPARK**:

- **proved free of run-time errors**: every index into memory, every
  multiplication in the indexing, every length;
- built on the **`blake2b` crate** for H0 and the variable-length hash H';
- its own **permutation P**, since Argon2's GB differs from BLAKE2b's by
  one multiplication (RFC 9106 §3.6);
- checked against **RFC 9106's test vectors** and against **libargon2** on
  random inputs;
- with the **PHC string** encoding used to store password hashes;
- wiping the memory it filled.

It is small, self-contained and useful far beyond any one application: a
proved password hash is what every Ada server that stores a password
needs.

---

## 2. The problem

### 2.1 Why Argon2id

RFC 9106 §4 names Argon2id the recommended variant and gives two parameter
sets:

- **FIRST RECOMMENDED**: t = 1, p = 4, m = 2^21 KiB (2 GiB), a 128-bit salt,
  a 256-bit tag;
- **SECOND RECOMMENDED**, for memory-constrained environments: t = 3,
  p = 4, m = 2^16 KiB (64 MiB), a 128-bit salt, a 256-bit tag.

### 2.2 What binding the C library costs

- a C dependency in the most sensitive path of a server (passwords);
- on some distributions the library is not in the base repositories, so
  it becomes a dependency carried by hand **[to verify per distribution]**;
- the guarantees a SPARK program has stop at the `Import`.

### 2.3 Why it suits SPARK

Argon2 is a fixed algorithm over arrays of 64-bit words: index
arithmetic, modular additions, rotations, one 32×32→64 multiplication,
and memory addressed by computed indices. Each of those is a run-time
check SPARK can discharge: the indices land inside memory, the
multiplication cannot overflow a 64-bit modular type, the lengths fit.
The structure is regular, so the proof should be tractable at a modest
level.

---

## 3. The algorithm, as RFC 9106 states it

### 3.1 Inputs (RFC §3.1)

| input | range |
|---|---|
| P, password | 0 .. 2^32-1 bytes |
| S, salt | up to 2^32-1 bytes; 16 bytes recommended, unique per password |
| p, parallelism (lanes) | 1 .. 2^24-1 |
| T, tag length | 4 .. 2^32-1 bytes |
| m, memory | 8·p .. 2^32-1 KiB |
| t, passes | 1 .. 2^32-1 |
| v, version | one byte, 0x13 |
| K, secret; X, associated data | optional, up to 2^32-1 bytes each |
| y, type | 0 Argon2d, 1 Argon2i, 2 Argon2id |

### 3.2 Steps

1. **H0** (§3.2): a 64-byte BLAKE2b hash over the parameters encoded
   little-endian, the lengths, and P, S, K, X.
2. **Memory** (§3.2, §3.4): m' = 4·p·⌊m / 4p⌋ blocks of 1024 bytes, a grid of
   p lanes by q = m'/p columns, cut into SL = 4 slices; a slice in one lane
   is a segment of q/4 blocks. Segments of one slice do not reference each
   other, so lanes can be filled in parallel.
3. **First blocks** of each lane from H' over H0, the column and the
   lane.
4. **Filling**: each block is G (previous block, reference block).
5. **Indexing for Argon2id** (§3.4.1.3, §3.4.2): in the **first pass,
   slices 0 and 1** use data-independent addressing (J1, J2 from a
   pseudo-random stream made of G calls, as Argon2i); **every other block**
   uses data-dependent addressing (J1, J2 from the previous block, as
   Argon2d). J2 mod p picks the reference lane (the current lane in the
   first slice of the first pass); J1 is mapped non-uniformly into the set
   of recently completed blocks (the RFC's Figures 12 and 13).
6. **G** (§3.5): R = X ⊕ Y; P applied to each row, then to each column, of
   the 8×8 matrix of 16-byte registers; the result XORed with R.
7. **P** (§3.6): the BLAKE2b round, with GB changed to
   `a = (a + b + 2 * trunc(a) * trunc(b)) mod 2^64` (trunc: the low
   32 bits), rotations 32, 24, 16, 63. **The multiplication is the only
   difference from BLAKE2b.**
8. **Finalisation** (§3.2 steps 7-8): the last column's blocks of every
   lane XORed into C; the tag is H'(C) of length T.

**H'** (§3.3): for T ≤ 64, one BLAKE2b call with output length T; for
longer T, a chain of 64-byte BLAKE2b outputs keeping 32 bytes of each, and
a last block sized to the remainder.

### 3.3 Test vectors (RFC §5.3)

Argon2id: 32 KiB, 3 passes, 4 lanes, a 32-byte tag; password 32 × 0x01,
salt 16 × 0x02, secret 8 × 0x03, associated data 12 × 0x04. The RFC gives
the pre-hashing digest, sample blocks after each pass, and the tag, so a
mismatch can be located to the step that diverges.

---

## 4. Design

### 4.1 Units

```
Argon2                  -- types, parameters, verdicts
Argon2.Blake2b_Bridge   -- H0 and H' over the blake2b crate
Argon2.Block            -- 1024-byte block as 128 Unsigned_64; XOR, copy
Argon2.Permutation      -- GB and P (SPARK, proved)
Argon2.Compress         -- G
Argon2.Index            -- J1/J2, the reference set, the non-uniform mapping
Argon2.Fill             -- passes, slices, segments
Argon2.Hash             -- the whole function
Argon2.PHC              -- $argon2id$v=19$m=..,t=..,p=..$salt$hash, encode and verify
Argon2.Wipe             -- memory cleared, not elidable
```

### 4.2 Memory: the caller's, or a managed wrapper

Argon2's memory is large (64 MiB to 2 GiB in the recommended sets), so it
cannot live on the stack. Two layers:

- **The SPARK core** takes the memory as a parameter:
  `Memory : in out Block_Array` with a precondition tying its length to m'.
  The core allocates nothing, which keeps it in SPARK without pointers.
- **A thin non-SPARK wrapper** allocates the array on the heap, calls the
  core, wipes the array and frees it. It is the API most users call.

### 4.3 The API

```ada
type Variant is (Argon2d, Argon2i, Argon2id);
type Parameters is record
   Kind       : Variant  := Argon2id;
   Passes     : Pass_Count;        -- t
   Memory_KiB : Memory_Size;       -- m
   Lanes      : Lane_Count;        -- p
   Tag_Length : Tag_Size;          -- T
end record;

type Verdict is (Ok, Parameters_Invalid, Below_Floor, Out_Of_Memory);

--  The RFC's whole range, for test vectors and interoperability.
procedure Hash_Raw (Params   : Parameters;
                    Password, Salt, Secret, Data : Byte_Array;
                    Memory   : in out Block_Array;
                    Tag      : out Byte_Array;
                    V        : out Verdict)
  with Pre  => Tag'Length = Params.Tag_Length
               and then Memory'Length = Blocks_For (Params),
       Post => Is_Wiped (Memory);

--  The everyday call: a PHC string, parameters checked against floors.
function Hash_Password (Password : Byte_Array; Params : Parameters := Recommended)
  return PHC_String;
function Verify (Password : Byte_Array; Stored : PHC_String) return Boolean;
--  Verify compares tags in constant time.
```

**Floors**: `Hash_Password` refuses parameters below a documented floor
(the SECOND RECOMMENDED set by default), so that a misconfiguration does
not silently store weak hashes. `Hash_Raw` accepts the RFC's whole range,
for test vectors and for interoperability with existing hashes.

### 4.4 Parallelism

Lanes of one slice are independent (§3.4), so filling them in parallel
is allowed and gives the same tag as filling them in turn. Version 1
**fills lanes sequentially** (correct for any p, simplest to prove).
Parallel filling comes later through tasks under the Jorvik profile, or
Ada 2022 parallel blocks **[to verify GNAT FSF's support]**, with a test
that the tag is identical.

### 4.5 Using the `blake2b` crate

H0 needs BLAKE2b with a 64-byte output; H' needs BLAKE2b with an output
length of 1 to 64 bytes, and an input that starts with a 4-byte length
prefix. The crate must offer **variable output length** and **streaming
input** **[to verify its API]**. If it does not, the bridge adds them,
upstream first, as a pull request to the crate.

The permutation P is **not** the crate's: GB's multiplication makes it a
different function (§3.6), so `Argon2.Permutation` implements it, reusing
the crate's rotation and round structure only as a reference for review.

---

## 5. Proof plan

| unit | properties | level |
|---|---|---|
| `Permutation` | no overflow (modular 64-bit), rotations within range | silver; functional against a ghost model of GB if tractable |
| `Index` | the 32×32→64 product and the shifts of the non-uniform mapping never overflow; the chosen column lies inside the reference set; the lane lies in 0 .. p-1 | silver plus postconditions |
| `Fill` | every block index computed lies inside `Memory`; the previous and the reference block exist when read; loops terminate | silver plus postconditions |
| `Blake2b_Bridge` | H' lengths: the chain covers T exactly; no index past a buffer | silver |
| `PHC` | the parser of a stored string is total: any string gives a result or a refusal | silver |
| `Wipe` | memory zero after the call | gold postcondition `Is_Wiped` |

No `pragma Annotate` justification without its reason written beside it,
and the proof replayed in CI at a fixed level.

---

## 6. Side channels

- **Argon2id leaks by design in its data-dependent half**: memory access
  depends on the password after the first half of the first pass. That is
  the variant's documented trade-off (data-independent first, so that a
  cache-timing observer learns little; data-dependent after, for
  resistance to GPU and trade-off attacks). The implementation adds no
  leak of its own.
- **The data-independent half** must be data-independent in the
  implementation too: no branch or index on the password there. Reviewed
  and measured (§7).
- **Tag comparison** in `Verify` is constant time.
- **Wiping**: the memory and every intermediate buffer (H0, the first
  blocks, C) are wiped through a volatile view the compiler cannot elide
  **[to verify the aspect's guarantee]**. A tag computation that fails
  half-way still wipes.

---

## 7. Testing

| test | how | oracle |
|---|---|---|
| RFC 9106 §5.3 vectors (Argon2id; also §5.1 and §5.2 for d and i) | exact, with the intermediate values the RFC prints | the RFC |
| differential on random inputs | password, salt, secret, data, t, m, p, T drawn at random within bounds; libargon2's C API (installed, 20190702) through a small C harness, which accepts K and X | libargon2 |
| the command line | `argon2 <salt> -id -t 3 -m 12 -p 4 -l 32 -r` (measured to work) for inputs without K and X | the `argon2` program |
| PHC strings | strings produced by the reference parsed and verified; ours accepted by the reference | the reference library |
| lane parallelism | sequential and parallel filling give the same tag | itself |
| fuzzing | the PHC parser and the parameter validation (companion project Harrow) | totality |
| timing | the data-independent phase measured with fixed and random passwords, in the style of dudect **[to verify the tool]** | statistics |
| performance | time and memory against libargon2 at the two recommended sets | numbers published; a target set after the first measurement, never before |

---

## 8. Performance

libargon2's reference has an optimised path using vector instructions.
A portable SPARK version will be slower at first. The aim, in order:

1. correct and proved;
2. measured against libargon2 at the SECOND RECOMMENDED set (64 MiB,
   t = 3, p = 4), with the ratio published;
3. made faster where the profile says so: `-O2`, unrolled loops in G,
   then parallel lanes, then vector code if GNAT can express it
   **[to verify: GCC vector extensions are not reachable from Ada
   directly; a C or assembly kernel would leave SPARK]**.

A password hash's slowness is partly its purpose: what matters is the
ratio to the attacker's hardware, not to a C implementation. A constant
factor is acceptable if it is measured and stated.

---

## 9. Roadmap

| milestone | content | done when |
|---|---|---|
| **M0** reading | RFC 9106 §3-§5 quoted into `docs/rfc9106-notes.md`; the PHC string format's specification read; the `blake2b` crate's API read | each [to verify] answered |
| **M1** H0 and H' | the bridge over `blake2b` | the RFC's pre-hashing digest reproduced exactly |
| **M2** G and P | `Permutation`, `Compress`, proved | a G test against values computed by libargon2's internals |
| **M3** indexing and filling | `Index`, `Fill`, sequential lanes | the RFC's intermediate blocks after pass 1 match; then the final tags of §5.1-§5.3 |
| **M4** proof | silver everywhere, the postconditions of §5 | gnatprove clean at the project's level, replayed in CI |
| **M5** API, PHC, wiping | §4.3, §6 | differential against libargon2 on 10^5 random inputs, PHC interoperability both ways |
| **M6** performance | measurement, then the first optimisations; parallel lanes | the ratio to libargon2 published |
| **M7** publication | crate, documentation, examples (a password store, a key derivation) | a stranger's CI reproduces the vectors and the proof |

---

## 10. Risks and open questions

| risk | answer |
|---|---|
| The `blake2b` crate lacks variable output or streaming | add it upstream; meanwhile a bridge |
| The proof of the indexing is harder than expected | the mapping is a handful of shifts and one product on 64-bit modular types; split into small lemmas |
| Too slow for the FIRST RECOMMENDED set (2 GiB) | documented; the SECOND set as the default; parallel lanes |
| Large allocations fail | `Out_Of_Memory` returned, never an exception escaping |
| A subtle divergence from the reference | the RFC's intermediate values locate it; the differential catches the rest |

Open questions for the maintainers:

1. Expose Argon2d and Argon2i publicly, or only Argon2id? (They share
   the code; exposing them helps interoperability, and invites misuse.)
2. Licence: Apache-2.0 (permissive, with a patent grant), or BSD-3-Clause
   to match the `blake2b` crate?
3. The default floor for `Hash_Password`: the SECOND RECOMMENDED set, or
   lower for embedded users with an explicit opt-in?

---

## 11. Repository layout

```
argon2id-spark/
  alire.toml                    -- depends on blake2b
  argon2id_spark.gpr
  src/
    argon2.ads
    argon2-blake2b_bridge.ad[sb]
    argon2-block.ad[sb]
    argon2-permutation.ad[sb]
    argon2-compress.ad[sb]
    argon2-index.ad[sb]
    argon2-fill.ad[sb]
    argon2-hash.ad[sb]
    argon2-phc.ad[sb]
    argon2-wipe.ad[sb]
    argon2-managed.ad[sb]       -- non-SPARK wrapper: allocate, call, wipe, free
  tests/
    rfc9106/                    -- vectors with intermediate values
    differential/               -- C harness over libargon2; CLI driver
    phc/
    timing/
  fuzz/                         -- Harrow targets
  proof/
  examples/
    password-store/  derive-key/
  docs/
    rfc9106-notes.md  proof.md  side-channels.md  performance.md
```

---

## Appendix A — measurements behind this document

21 October 2026, Linux x86-64:

```bash
alr search --crates | grep -i argon              # nothing
alr show blake2b                                  # 0.1.1, BSD-3-Clause, tags spark/proof
dpkg -l libargon2-1 libargon2-dev                 # 0~20190702+dfsg-4build1
ls /usr/include/argon2.h                          # present
echo -n password | argon2 somesalt16bytes! -id -t 3 -m 12 -p 4 -l 32 -r
#   c750b3329755075a658f5beab9c839b9bd9278b7034df66b4fece272a265ff58
gnatprove --version                               # FSF 16.1.0
```

RFC 9106 read from <https://www.rfc-editor.org/rfc/rfc9106.txt>.
