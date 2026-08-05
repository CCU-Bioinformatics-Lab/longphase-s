// Embeds the somatic indel CNN weights into the binary, so longphase-s needs no
// companion model file at runtime.
//
// The weights live in src/somatic_haplotag/{del,ins}_model.bin, exported from
// the training checkpoints:
//   del_model.bin  <- output_multi_FA_event_thr_0.5_del_only/best_model.pt
//   ins_model.bin  <- output_multi_FA_event_thr_0.5_ins_only_0.12budget/best_model.pt
// To ship a retrained model, replace the relevant file and rebuild; nothing in
// this file changes.
//
// The .incbin paths are resolved by the assembler relative to the compiler's
// working directory, which the Makefile always sets to the repository root.

asm(".section .rodata\n"
    ".global SOMATIC_DEL_MODEL_WEIGHTS\n"
    ".balign 4\n"
    "SOMATIC_DEL_MODEL_WEIGHTS:\n"
    ".incbin \"src/somatic_haplotag/del_model.bin\"\n"
    ".global SOMATIC_DEL_MODEL_WEIGHTS_END\n"
    "SOMATIC_DEL_MODEL_WEIGHTS_END:\n"

    ".global SOMATIC_INS_MODEL_WEIGHTS\n"
    ".balign 4\n"
    "SOMATIC_INS_MODEL_WEIGHTS:\n"
    ".incbin \"src/somatic_haplotag/ins_model.bin\"\n"
    ".global SOMATIC_INS_MODEL_WEIGHTS_END\n"
    "SOMATIC_INS_MODEL_WEIGHTS_END:\n"
    ".previous\n");
