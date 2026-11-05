const repositoryUrl = "https://github.com/StatFunGen/P8149_2026_vignettes";
const sourceFiles = {
  "p8149-lecture-2-em-basics": "P8149_Lecture_2_EM_Basics.ipynb",
  "p8149-lecture-2-maf-mle-em": "P8149_Lecture_2_MAF_MLE_EM.ipynb",
  "p8149-lecture-2-hmm-basics": "P8149_Lecture_2_HMM_Basics.ipynb",
  "p8149-lecture-2-li-stephens-hmm": "P8149_Lecture_2_Li_Stephens_HMM.ipynb",
  "p8149-lecture-3-local-ibd": "P8149_Lecture_3_Local_IBD.ipynb",
  "p8149-lecture-3-genome-wide-kinship": "P8149_Lecture_3_Genome_Wide_Kinship.ipynb",
};

function render({ el }) {
  const style = document.createElement("style");
  style.textContent = [
    ":host { display: inline-flex; align-items: center; }",
    ".source-link { display: inline-flex; align-items: center; padding: 0.42rem 0.82rem; border: 1px solid #c8cfd8; border-radius: 0.42rem; color: inherit; background: transparent; font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; font-size: 0.9rem; font-weight: 600; line-height: 1.15; text-decoration: none; white-space: nowrap; }",
    ".source-link:hover { border-color: #8b96a3; background: rgba(127, 127, 127, 0.08); }",
  ].join("\n");

  const link = document.createElement("a");
  link.className = "source-link";
  link.href = repositoryUrl;
  link.target = "_blank";
  link.rel = "noopener noreferrer";
  link.textContent = "GitHub";
  link.setAttribute("aria-label", "View this page on GitHub");

  el.append(style, link);

  let currentSource = "";
  const updateLink = () => {
    const editLink = document.querySelector('a[aria-label="Edit This Page"]');
    const slug = window.location.pathname.replace(/\/$/, "").split("/").pop();
    const sourceFile = sourceFiles[slug] || "index.md";
    const source = editLink?.href?.replace("/edit/", "/blob/") ||
      repositoryUrl + "/blob/main/" + sourceFile;
    if (source !== currentSource) {
      link.href = source;
      currentSource = source;
    }
  };

  updateLink();
  const timer = window.setInterval(updateLink, 500);

  return () => {
    window.clearInterval(timer);
    link.remove();
    style.remove();
  };
}

export default { render };
