import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG
import json
import base64
from io import BytesIO

st.set_page_config(layout="wide")

# Helper to draw subgraph with atom indices
def draw_with_atom_indices(mol, atom_indices=None):
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.drawOptions().addAtomIndices = True
    if atom_indices:
        drawer.DrawMolecule(mol, highlightAtoms=atom_indices)
    else:
        drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')
    return svg

# Page title
st.title("🧬 Molecular Annotation UI")

# Sidebar navigation
mode = st.sidebar.radio("Choose annotation mode:", ["Training Set", "Test Set"])

# Shared: Input SMILES
st.sidebar.markdown("### Input Molecules")
train_smiles_list = st.sidebar.text_area("Training Molecules (SMILES, one per line):",
                                         "CCO\nCCN\nCNC")
test_smiles_list = st.sidebar.text_area("Test Molecules (SMILES, one per line):",
                                        "CCO\nCCCl")

train_molecules = [s.strip() for s in train_smiles_list.split('\n') if s.strip()]
test_molecules = [s.strip() for s in test_smiles_list.split('\n') if s.strip()]

if mode == "Training Set":
    st.header("🔧 Training Molecule Annotation")

    selected = st.selectbox("Select molecule to annotate:", train_molecules)
    mol = Chem.MolFromSmiles(selected)
    if mol:
        st.subheader("Molecule Viewer (Atom Indices)")
        svg = draw_with_atom_indices(mol)
        st.components.v1.html(svg, height=300)

        st.subheader("Subgraph Annotation")
        subgraph_count = st.number_input("How many subgraphs to annotate?", min_value=1, max_value=10, value=2)

        annotations = []
        for i in range(subgraph_count):
            st.markdown(f"**Subgraph {i+1}**")
            col1, col2 = st.columns([1, 2])
            with col1:
                st.markdown("Original Molecule")
                st.components.v1.html(svg, height=300)
            with col2:
                atom_indices_str = st.text_input(f"Atom indices for subgraph {i+1} (comma-separated)", key=f"idx_{i}")
                weight = st.slider(f"Importance weight for Subgraph {i+1}", 0.0, 1.0, 0.5, 0.05, key=f"w_{i}")
                if atom_indices_str:
                    try:
                        atom_indices = [int(x) for x in atom_indices_str.split(',') if x.strip().isdigit()]
                        sub_svg = draw_with_atom_indices(mol, atom_indices)
                        st.components.v1.html(sub_svg, height=300)
                        atom_indices = [f"atom_{idx}" for idx in atom_indices] # Add atom in front of idx to avoid confusion
                        annotations.append({"subgraph": atom_indices, "weight": weight})
                    except Exception as e:
                        st.warning(f"Invalid atom indices: {e}")

        if st.button("💾 Save Annotations for This Molecule"):
            st.success(f"Annotations saved for molecule: {selected}")
            annotation_data = {"molecule": selected, "annotations": annotations}
            st.json(annotation_data)
            json_str = json.dumps(annotation_data, indent=2)
            b64 = base64.b64encode(json_str.encode()).decode()
            href = f'<a href="data:file/json;base64,{b64}" download="annotation_train_{selected}.json">📥 Download annotation file</a>'
            st.markdown(href, unsafe_allow_html=True)

elif mode == "Test Set":
    st.header("🔍 Test Molecule Pairwise Ranking")

    test_selected = st.selectbox("Select test molecule:", test_molecules)
    test_mol = Chem.MolFromSmiles(test_selected)

    if test_mol:
        st.subheader("Test Molecule Viewer (Atom Indices)")
        svg = draw_with_atom_indices(test_mol)
        st.components.v1.html(svg, height=300)

        st.subheader("Pairwise Ranking Against Training Set")
        pairwise_annotations = []
        for ref_smiles in train_molecules:
            ref_mol = Chem.MolFromSmiles(ref_smiles)
            col1, col2, col3 = st.columns([1, 1, 2])

            with col1:
                st.markdown("**Test**")
                st.image(Draw.MolToImage(test_mol, size=(300, 300)))
            with col2:
                st.markdown("**Reference**")
                st.image(Draw.MolToImage(ref_mol, size=(300, 300)))
            with col3:
                choice = st.radio(f"Which has higher property value? (Test vs {ref_smiles})",
                                  ["Test > Reference", "Test < Reference", "Equal / Uncertain"],
                                  key=f"rank_{ref_smiles}")
                pairwise_annotations.append({"reference": ref_smiles, "comparison": choice})

        st.subheader("Subgraph Importance for Test Molecule")
        subgraph_count_test = st.number_input("How many subgraphs to annotate (test)?", min_value=1, max_value=10, value=2, key="test_subgraphs")

        test_annotations = []
        for i in range(subgraph_count_test):
            st.markdown(f"**Subgraph {i+1}**")
            col1, col2 = st.columns([1, 2])
            with col1:
                st.markdown("Original Molecule")
                st.components.v1.html(svg, height=300)
            with col2:
                atom_indices_str = st.text_input(f"Atom indices for subgraph {i+1} (comma-separated)", key=f"t_idx_{i}")
                weight = st.slider(f"Importance weight for Subgraph {i+1}", 0.0, 1.0, 0.5, 0.05, key=f"t_w_{i}")
                if atom_indices_str:
                    try:
                        atom_indices = [int(x) for x in atom_indices_str.split(',') if x.strip().isdigit()]
                        sub_svg = draw_with_atom_indices(test_mol, atom_indices)
                        st.components.v1.html(sub_svg, height=300)
                        atom_indices = [f"atom_{idx}" for idx in atom_indices] # Add atom in front of idx to avoid confusion
                        test_annotations.append({"subgraph": atom_indices, "weight": weight})
                    except Exception as e:
                        st.warning(f"Invalid atom indices: {e}")

        if st.button("💾 Save Annotations for This Test Molecule"):
            st.success(f"Annotations saved for test molecule: {test_selected}")
            annotation_data = {"molecule": test_selected, "pairwise": pairwise_annotations, "annotations": test_annotations}
            st.json(annotation_data)
            json_str = json.dumps(annotation_data, indent=2)
            b64 = base64.b64encode(json_str.encode()).decode()
            href = f'<a href="data:file/json;base64,{b64}" download="annotation_test_{test_selected}.json">📥 Download annotation file</a>'
            st.markdown(href, unsafe_allow_html=True)
