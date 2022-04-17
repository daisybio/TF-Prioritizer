import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactor} from "../../types/types";

@Component({
  selector: 'single-tf-accordion',
  templateUrl: './single-tf-accordion.component.html',
  styleUrls: ['./single-tf-accordion.component.css']
})
export class SingleTfAccordionComponent implements OnInit {
  @Input()
  // @ts-ignore
  transcriptionFactor: TranscriptionFactor;

  visible: boolean = false;

  constructor() { }

  ngOnInit(): void {
  }

  toggleVisibility() {
    this.visible = !this.visible;
  }
}
