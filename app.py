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

img_width = 800
img_height = 500

# Helper to draw subgraph with atom indices
def draw_with_atom_indices(mol, atom_indices=None):
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(img_width, img_height)
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
                                        "C#CC1=C2C(=O)[Si](C(C)=O)(C(=O)OC)C(=O)C2=C(C#Cc2ccc(-c3nnc(-c4cccs4)n3C(=O)OC)s2)C1(C(C)=O)C(=O)OC\n\
                                        COc1c2scc(C(F)F)c2c(OCF)c2sc(-c3ccc(-c4sc(-c5cccs5)c5cc(C(=O)N(C)C)c(C(O)N(C)C)cc45)s3)c(C(F)(F)F)c12\n\
                                        CC(C)(C)c1cc(CC2(C=C(C#N)C#N)C(=O)c3csc(-c4ccc(C5=C6C(=O)N(C#N)C(c7cccs7)=C6C(=O)N5C=C(C#N)C#N)s4)c3C2=O)cc(C(C)(C)C)c1\n\
                                        CNc1ccc2cc(NCc3ccc(-c4cc5c(s4)C4=[SH]C(c6cccn6C)=CC4=C5C(C)C)n3C)ccc2c1\n\
                                        Cc1c(-c2nnc(-c3cc4sccc4s3)o2)sc(-c2nnc(-c3cc4sc(-c5ccc6c(c5)n(C(C)C)c5c7ccccc7n(C(C)C)c65)cc4s3)o2)c1C\n\
                                        CO[Si]1(C)C(=O)c2c(c(-c3sc(-c4ccc(-c5nnc(-c6cccs6)n5C)s4)c4c3OCCO4)c3c(c2-c2scc4c2OCCO4)C(=O)N(C(C)=O)C3=O)C1=O")

train_molecules = [s.strip() for s in train_smiles_list.split('\n') if s.strip()]
test_molecules = [s.strip() for s in test_smiles_list.split('\n') if s.strip()]

if mode == "Training Set":
    st.header("🔧 Training Molecule Annotation")

    selected = st.selectbox("Select molecule to annotate:", train_molecules)
    mol = Chem.MolFromSmiles(selected)
    if mol:
        st.subheader("Molecule Viewer (Atom Indices)")
        svg = draw_with_atom_indices(mol)
        st.components.v1.html(svg, height=img_height, scrolling=True)

        st.subheader("Subgraph Annotation")
        subgraph_count = st.number_input("How many subgraphs to annotate?", min_value=1, max_value=10, value=2)

        annotations = []
        for i in range(subgraph_count):
            st.markdown("---")
            st.markdown(f"**Subgraph {i+1}**")
            
            atom_indices_str = st.text_input(f"Atom indices for subgraph {i+1} (comma-separated)", key=f"idx_{i}")
            weight = st.slider(f"Importance weight for Subgraph {i+1}", 0.0, 1.0, 0.5, 0.05, key=f"w_{i}")
            try:
                atom_indices = [int(x) for x in atom_indices_str.split(',') if x.strip().isdigit()]
                sub_svg = draw_with_atom_indices(mol, atom_indices)
                st.components.v1.html(sub_svg, height=img_height, scrolling=True)
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

    st.sidebar.markdown("### Properties of Interest")
    properties_input = st.sidebar.text_input("Enter properties (comma-separated):", "Activity, Solubility")
    properties = [p.strip() for p in properties_input.split(',') if p.strip()]

    test_selected = st.selectbox("Select test molecule:", test_molecules)
    test_mol = Chem.MolFromSmiles(test_selected)

    if test_mol:
        st.subheader("Test Molecule Viewer (Atom Indices)")
        svg = draw_with_atom_indices(test_mol)
        st.components.v1.html(svg, height=img_height, scrolling=True)

        st.subheader("Pairwise Ranking Against Training Set")
        pairwise_annotations = []
        for ref_smiles in train_molecules:
            ref_mol = Chem.MolFromSmiles(ref_smiles)
            col1, col2, col3 = st.columns([1, 1, 2])

            with col1:
                st.markdown("**Test**")
                st.image(Draw.MolToImage(test_mol, size=(img_width, img_height)))
            with col2:
                st.markdown("**Reference**")
                st.image(Draw.MolToImage(ref_mol, size=(img_width, img_height)))
            with col3:
                prop_choices = {}
                for prop in properties:
                    prop_choice = st.radio(f"{prop}: Test vs {ref_smiles}",
                                           ["Test > Reference", "Test < Reference", "Equal / Uncertain"],
                                           key=f"rank_{ref_smiles}_{prop}")
                    prop_choices[prop] = prop_choice
                pairwise_annotations.append({"reference": ref_smiles, "comparisons": prop_choices})
            st.markdown("---")


        st.subheader("Subgraph Importance for Test Molecule")
        subgraph_count_test = st.number_input("How many subgraphs to annotate (test)?", min_value=1, max_value=10, value=2, key="test_subgraphs")

        test_annotations = []
        for i in range(subgraph_count_test):
            st.markdown("---")
            st.markdown(f"**Subgraph {i+1}**")
            
            atom_indices_str = st.text_input(f"Atom indices for subgraph {i+1} (comma-separated)", key=f"t_idx_{i}")
            weights = {}
            for prop in properties:
                weights[prop] = st.slider(f"Importance for {prop} (Subgraph {i+1})", 0.0, 1.0, 0.5, 0.05, key=f"t_w_{i}_{prop}")
            try:
                atom_indices = [int(x) for x in atom_indices_str.split(',') if x.strip().isdigit()]
                sub_svg = draw_with_atom_indices(test_mol, atom_indices)
                st.components.v1.html(sub_svg, height=img_height, scrolling=True)
                atom_indices = [f"atom_{idx}" for idx in atom_indices] # Add atom in front of idx to avoid confusion
                test_annotations.append({"subgraph": atom_indices, "weights": weights})
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
