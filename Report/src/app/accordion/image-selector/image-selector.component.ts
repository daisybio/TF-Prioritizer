import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'image-selector',
  templateUrl: './image-selector.component.html',
  styleUrls: ['./image-selector.component.css']
})
export class ImageSelectorComponent implements OnInit {
  @Input()
    // @ts-ignore
  data: {};

  @Input()
  visible: boolean = false;

  levels: Set<string>[] = [];
  activeLevels: string[] = [];

  hasData: boolean = false;
  hasPlots: boolean = false;

  imageSource: string = "";
  dataPrefix: string = "assets/input/data"

  constructor() {
  }

  ngOnInit(): void {
    this.collectLevels(this.data, 0);

    this.hasData = this.levels[this.levels.length - 1].has("data");
    this.hasPlots = this.levels[this.levels.length - 1].has("plot");

    this.levels.pop();
    this.activeLevels.pop();

    this.updateImage();
  }

  updateImage() {
    let subMap = this.data;

    for (let level of this.activeLevels) {
      // @ts-ignore
      subMap = subMap[level];
    }

    // @ts-ignore
    this.imageSource = this.dataPrefix + subMap["plot"];
  }

  downloadData() {
    let subMap = this.data;

    for (let level of this.activeLevels) {
      // @ts-ignore
      subMap = subMap[level];
    }

    const link = document.createElement("a");
    link.setAttribute("target", "_blank");
    // @ts-ignore
    let file: string = subMap["data"];
    link.setAttribute("href", this.dataPrefix + file);
    link.setAttribute("download", "test.csv");
    document.body.appendChild(link);
    link.click();
    link.remove();
  }

  collectLevels(map: {}, level: number) {
    if (this.levels.length <= level) {
      this.levels.push(new Set());
    }

    for (let key of Object.keys(map)) {
      this.levels[level].add(key);
      if (this.activeLevels.length <= level) {
        this.activeLevels.push(key);
      }

      // @ts-ignore
      let sub = map[key];
      if (!(sub instanceof String || typeof sub == 'string')) {
        this.collectLevels(sub, level + 1);
      }
    }
  }
}
