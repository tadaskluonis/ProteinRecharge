# ProteinRecharge Tolerance Parameter Changes

## Summary
Replaced single `tolerance` parameter with two separate tolerance parameters for finer control over charge optimization:
- `batch_charge_tolerance`: Controls when to stop bias adjustments
- `charge_tolerance`: Controls acceptable charge range for final sequence selection

## Changes Made

### Config Dataclass (recharge.py)
**Before:**
```python
desired_chain_charge: float = 0.0
tolerance: float = 1.0
initial_bias: float = 0.0
```

**After:**
```python
desired_chain_charge: float = 0.0
batch_charge_tolerance: float = 1.0  # Tolerance for mean batch charge to stop bias adjustment
charge_tolerance: float = 1.0        # Tolerance for individual sequence charge in final pick
initial_bias: float = 0.0
```

### Run Loop Logic (run_recharge_pipeline)
The iteration loop now implements a two-stage process:

1. **Bias Adjustment Phase:**
   - Continue adjusting bias until mean batch charge is within `batch_charge_tolerance` of desired charge

2. **Sequence Selection Phase:**
   - Once batch tolerance is reached, filter sequences by `charge_tolerance`
   - If sequences meet criteria: pick one with best `overall_confidence`
   - If no sequences meet criteria: continue bias adjustment with new batch

**Old Logic:**
```python
if abs(mean_q - desired_charge) <= tolerance:
    print("Tolerance reached.")
    break
# Then pick closest charge from final batch
```

**New Logic:**
```python
if abs(mean_q - desired_charge) <= batch_charge_tolerance:
    print(f"Batch charge tolerance reached. Looking for sequence with charge within ±{charge_tolerance}...")
    candidates = [r for r in current_results if abs(r['charge'] - desired_charge) <= charge_tolerance]
    
    if candidates:
        best_pick = max(candidates, key=lambda x: x.get('overall_confidence', 0))
        break
    else:
        print(f"No sequences found within charge tolerance ±{charge_tolerance}. Continuing bias adjustment...")
```

### YAML Config Update
**Before:**
```yaml
desired_chain_charge: 6.0
tolerance: 4.0
```

**After:**
```yaml
desired_chain_charge: 6.0
batch_charge_tolerance: 4.0  # Tolerance for mean batch charge to stop bias adjustment
charge_tolerance: 1.0        # Tolerance for individual sequence charge in final pick
```

## Parameter Guidance

- **batch_charge_tolerance**: Should be relatively loose to ensure the algorithm reaches a point where it's testing sequences near the target charge. Typical range: 0.5-2.0 per charge unit
- **charge_tolerance**: Should be tighter to ensure selected sequences are within your desired range. Typical range: 0.5-1.0 per charge unit

Example:
```yaml
desired_chain_charge: 5.0
batch_charge_tolerance: 1.5  # Stop adjusting once mean batch is within ±1.5 of target
charge_tolerance: 0.5        # Pick sequence with charge within ±0.5 of target
```

## Benefits

1. **More Control**: Fine-tune convergence (batch tolerance) separately from final sequence acceptance (charge tolerance)
2. **Better Selectivity**: When sequences meet the tighter `charge_tolerance`, pick the best confidence one rather than just the closest
3. **Continued Search**: If no sequences meet criteria after reaching batch tolerance, automatically continue searching rather than stopping
4. **Explicit Intent**: Parameter names clearly indicate their purpose in the optimization workflow
