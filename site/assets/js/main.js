(() => {
  const header = document.querySelector("[data-header]");
  if (header) {
    const onScroll = () => {
      header.classList.toggle("is-scrolled", window.scrollY > 4);
    };
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
  }

  const reduce = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  document.querySelectorAll("[data-copy]").forEach((btn) => {
    btn.addEventListener("click", async () => {
      const sel = btn.getAttribute("data-copy");
      const node = sel ? document.querySelector(sel) : null;
      if (!node) return;
      const text = node.textContent || "";
      try {
        await navigator.clipboard.writeText(text);
      } catch {
        const range = document.createRange();
        range.selectNodeContents(node);
        const pick = window.getSelection();
        pick.removeAllRanges();
        pick.addRange(range);
      }
      const status = document.querySelector("[data-copy-status]");
      if (status) {
        status.textContent = "Copied.";
        status.hidden = false;
        window.setTimeout(() => {
          if (!reduce) status.hidden = true;
        }, 1600);
      }
    });
  });
})();
