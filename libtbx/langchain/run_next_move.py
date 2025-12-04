from __future__ import division
import asyncio
from libtbx.langchain import langchain_tools as lct
import sys

def run():
    # --- 1. Setup ---
    try:
        print("Initializing Agent...")
        # Use the PRO model for the "Brain" of the operation
        llm, embeddings = lct.get_llm_and_embeddings( provider='google',
           llm_model_name='gemini-2.5-pro'
         )

    except Exception as e:
        print(f"Setup failed: {e}")
        return

    # --- 2. Load History (Mocked for this example) ---
    # In production, you would read these from your 'log_info' objects or files
    history = [
        {
            "summary": "Run 1: phenix.refine on model.pdb and data.mtz. R-work: 0.48, R-free: 0.52. Geometry is poor.",
            "analysis": "The R-factors are very high. The model likely needs significant conformational changes. Rigid body refinement was not sufficient."
        },
        {
            "summary": "Run 2: phenix.refine with rigid_body strategy. R-work: 0.45, R-free: 0.49. Slight improvement.",
            "analysis": "Better, but still stuck in a local minimum. The structure might need simulated annealing to overcome energy barriers."
        }
    ]

    # --- 3. Run Super-Analysis ---
    print("\n--- Starting Super-Analysis ---")
    result = asyncio.run(lct.generate_next_move(
        run_history=history,
        llm=llm,
        embeddings=embeddings,
        db_dir='/net/cci-gpu-01/raid1/scratch1/terwill/build_dir/modules/phenix/phenix/phenix_ai/docs_db_google',
    ))

    # --- 4. Output Results ---
    if result.error:
        print(f"\nERROR: {result.error}")
    else:
        print("\n" + "="*40)
        print("RECOMMENDED ACTION")
        print("="*40)
        print(f"**Explanation:**\n{result.explanation}\n")
        print(f"**Strategy:**\n{result.strategy}\n")
        print("="*40)
        print("GENERATED COMMAND:")
        print(f"{result.command}")
        print("="*40)

if __name__ == "__main__":
    run()
