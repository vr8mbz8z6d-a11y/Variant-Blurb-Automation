"""
Streamlit web front end for the variant blurb pipeline.

This does NOT reimplement any pipeline logic -- it's a thin UI wrapper
around the exact same functions the command-line script
(Variant_Blurb_Automation.py) uses: parse_combined_variant_string,
resolve_coordinates, annotate_and_render, normalize_clinical_classification.
If a bug gets fixed in the underlying pipeline modules, this app picks up
the fix automatically, with nothing here to keep in sync.

DEPLOYMENT (Streamlit Community Cloud, free tier):
  1. Push this whole project folder to a GitHub repo (private repo is
     fine -- Streamlit Cloud can deploy from private repos once you
     connect your GitHub account).
  2. Go to https://share.streamlit.io, sign in with GitHub, click
     "New app", point it at this repo and this file (app.py).
  3. Before (or after) deploying, set the password secret: in the
     Streamlit Cloud dashboard for this app, go to Settings -> Secrets,
     and paste:
         app_password = "choose-a-password-here"
     This is NOT committed to your GitHub repo -- it's stored separately
     by Streamlit Cloud, which is exactly why it's safe to put a real
     password there (unlike hardcoding it in this file).
  4. Share the resulting *.streamlit.app URL with your two colleagues,
     along with the password out-of-band (e.g. a text/Slack message,
     not written anywhere in the repo itself). From their side, this
     behaves exactly like clicking a link and using a normal web page --
     no Python, no installation, nothing to run locally.

PRIVACY NOTE, stated plainly rather than implied: this password check is
a basic access gate, not strong security -- it stops casual/accidental
access by someone who finds the link, but a determined party could
still get in. Given the tool only ever handles variant-level HGVS
strings (never patient names/identifiers), this is treated as an
appropriately-scoped protection for the actual sensitivity of what's
being entered here -- not a claim that this meets any formal compliance
bar. If that changes (e.g. if patient-identifying info were ever
entered), a stronger access-control approach would be needed.
"""
import io
import contextlib
import streamlit as st

from Variant_Blurb_Automation import annotate_and_render, HGVSResolutionError
from ensembl_hgvs_source import resolve_coordinates
from variant_input_parser import parse_combined_variant_string, VariantStringParseError
from templates import normalize_clinical_classification

st.set_page_config(page_title="Variant Blurb Tool", page_icon="🧬", layout="centered")


def _check_password() -> bool:
    """
    Simple password gate. Reads the expected password from Streamlit's
    secrets store (st.secrets["app_password"]) -- never hardcoded here,
    so it's safe for this file to live in a shared/private repo.
    """
    if "authenticated" not in st.session_state:
        st.session_state.authenticated = False

    if st.session_state.authenticated:
        return True

    st.title("🧬 Variant Blurb Tool")
    password = st.text_input("Password", type="password")
    if st.button("Enter"):
        expected = st.secrets.get("app_password")
        if expected is None:
            st.error(
                "No app_password is configured in Streamlit secrets yet. "
                "See the deployment instructions in app.py's docstring."
            )
        elif password == expected:
            st.session_state.authenticated = True
            st.rerun()
        else:
            st.error("Incorrect password.")
    return False


def main():
    if not _check_password():
        return

    st.title("🧬 Variant Blurb Tool")
    st.caption(
        "Enter a variant as a single combined string, e.g.: "
        "`NM_000326.5(RLBP1):c.753C>A; p.(Tyr251*)`"
    )
    st.caption(
        "Transcript accession is required. Gene symbol and p. notation are "
        "optional -- splice variants typically have no p. notation."
    )

    raw_input_str = st.text_input("Variant", placeholder="NM_000326.5(RLBP1):c.753C>A; p.(Tyr251*)")
    genome_build = st.selectbox("Genome build", ["hg38", "hg19"], index=0)
    classification_input = st.text_input(
        "Variant classification (optional)",
        placeholder="e.g. Pathogenic, Likely benign, VUS -- leave blank to use ClinVar's own call",
        help="If given, this overrides ClinVar's aggregate classification in the closing "
             "summary sentence. Leave blank to let the pipeline use ClinVar's own call "
             "automatically.",
    )

    with st.expander("Advanced"):
        show_debug = st.checkbox(
            "Show debug output (raw API requests/responses)", value=False,
            help="Useful for diagnosing an unexpected result -- shows exactly what each "
                 "external database returned.",
        )

    if st.button("Generate blurb", type="primary"):
        if not raw_input_str.strip():
            st.warning("Enter a variant first.")
            return

        try:
            parsed = parse_combined_variant_string(raw_input_str)
        except VariantStringParseError as e:
            st.error(f"Could not parse that input: {e}")
            return

        clinical_classification = normalize_clinical_classification(classification_input.strip() or None)
        if classification_input.strip() and clinical_classification is None:
            st.warning(
                f"Unrecognized classification {classification_input.strip()!r} -- "
                f"proceeding with ClinVar's own classification instead."
            )

        st.write(
            f"**Parsed** -- transcript: `{parsed.transcript}`, "
            f"gene: `{parsed.gene or '(not given)'}`, "
            f"hgvs_c: `{parsed.hgvs_c}`, "
            f"hgvs_p: `{parsed.hgvs_p or '(not given)'}`"
        )

        debug_capture = io.StringIO()
        try:
            with st.spinner("Querying Ensembl, gnomAD, ClinVar, and related databases..."):
                with contextlib.redirect_stdout(debug_capture):
                    # Resolve coordinates FIRST and separately, exactly as
                    # the interactive script does -- see
                    # Variant_Blurb_Automation.py's __main__ block for the
                    # full reasoning (avoids ever running downstream
                    # lookups against placeholder coordinates if
                    # resolution fails).
                    resolved = resolve_coordinates(
                        parsed.hgvs_c_with_transcript, genome_build=genome_build, debug=show_debug
                    )
                    record, blurb = annotate_and_render(
                        chrom=resolved.chrom, pos=resolved.pos, ref=resolved.ref, alt=resolved.alt,
                        hgvs_c=parsed.hgvs_c, hgvs_p=parsed.hgvs_p, gene=parsed.gene,
                        genome_build=genome_build,
                        hgvs_c_with_transcript=parsed.hgvs_c_with_transcript,
                        clinical_classification=clinical_classification,
                        debug=show_debug,
                    )
        except HGVSResolutionError as e:
            st.error(
                f"Could not resolve genomic coordinates for this variant:\n\n{e}\n\n"
                f"Common causes: malformed HGVS notation (e.g. missing the '>alt' half "
                f"of a substitution), a transient network issue (try again), or a "
                f"transcript accession that doesn't match the gene. Nothing downstream "
                f"was queried, so no confusing blurb was built from wrong coordinates."
            )
            if show_debug and debug_capture.getvalue():
                with st.expander("Debug output"):
                    st.code(debug_capture.getvalue())
            return
        except Exception as e:
            # Catch-all so a shared, multi-person tool never shows a raw
            # Python traceback to a colleague -- always something
            # readable, with the debug log available if they turned it on.
            st.error(f"Something went wrong: {e}")
            if show_debug and debug_capture.getvalue():
                with st.expander("Debug output"):
                    st.code(debug_capture.getvalue())
            return

        st.subheader("Result")
        st.write(f"**Variant type:** {record.variant_type}")
        st.write(f"**Resolved transcript:** {record.hgvs_transcript}")
        st.markdown("**Blurb:**")
        if blurb:
            st.info(blurb)
        else:
            st.warning("No facts resolved -- check network access / coordinates.")

        if show_debug and debug_capture.getvalue():
            with st.expander("Debug output", expanded=False):
                st.code(debug_capture.getvalue())


if __name__ == "__main__":
    main()
