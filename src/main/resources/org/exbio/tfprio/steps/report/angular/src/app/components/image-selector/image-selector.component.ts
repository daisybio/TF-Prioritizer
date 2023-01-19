import {Component, Input} from '@angular/core';
import {Subject} from "rxjs";

@Component({
  selector: 'app-image-selector',
  templateUrl: './image-selector.component.html',
  styleUrls: ['./image-selector.component.scss']

})
export class ImageSelectorComponent {
  @Input()
  data: {} = {};
  activeData: Subject<string> = new Subject<string>();
}

